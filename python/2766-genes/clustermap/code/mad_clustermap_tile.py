#!/usr/bin/env python
# coding: utf-8

# %% cml input
import sys, os

# tile_num = int(sys.argv[1])
# tile_num = 315
# tile = f"Position{tile_num:03}"

tile = sys.argv[1]
print(tile)

from ClusterMap.clustermap import *
from anndata import AnnData
import tifffile
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import scipy.io
from sklearn import preprocessing
import timeit
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy import ndimage
import copy, math

start = timeit.default_timer()

# %%
### set file folder
ppath = '/stanley/WangLab/Data/Processed/2022-01-03-Hu-AD_8m_Rep/output/max/'
dpath = '/stanley/WangLab/Data/Processed/2022-01-03-Hu-AD_8m_Rep/output/max/pi/'
opath = '/stanley/WangLab/Data/Processed/2022-01-03-Hu-AD_8m_Rep/output/clustermap/'

if not os.path.exists(opath):
    os.mkdir(opath)

### set parameters
xy_radius = 50 # for knn radius (only)?
z_radius = 10
pct_filter = 0.1
dapi_grid_interval = 5
min_spot_per_cell = 10
cell_num_threshold = 0.02
window_size = 512

# %%
### read dapi: col, row, z
dapi = tifffile.imread(os.path.join(dpath, f"{tile}.tif"))
dapi = np.transpose(dapi, (1,2,0))

### read spots
mat = scipy.io.loadmat(os.path.join(ppath, tile, 'goodPoints_max3d.mat'))
spots = pd.DataFrame(mat['tile_goodSpots'])
spots.columns=['spot_location_1','spot_location_2','spot_location_3']

### convert gene_name to gene identity
gene_name = list(x[0][0] for x in mat['tile_goodReads'])
le = preprocessing.LabelEncoder()
le.fit(gene_name)
gene = le.transform(gene_name) + 1 #minimal value of gene=1

### save genelist
genes = pd.DataFrame(le.classes_)
genes_opath = os.path.join(opath, tile)
if not os.path.exists(genes_opath):
    os.mkdir(genes_opath)
genes.to_csv(os.path.join(genes_opath, 'genelist.csv'), header=False, index=False)
spots['gene'] = gene
spots['gene'] = spots['gene'].astype('int')

### instantiate model
num_gene = np.max(spots['gene'])
gene_list = np.arange(1, num_gene+1)
num_dims = len(dapi.shape)
model = ClusterMap(spots=spots, dapi=dapi, gene_list=gene_list, num_dims=num_dims, # gauss_blur=True, sigma=8,
               xy_radius=xy_radius, z_radius=z_radius, fast_preprocess=False)

model.preprocess(dapi_grid_interval=dapi_grid_interval, pct_filter=pct_filter)

# Since tiles are large, split data into smaller tiles and re-stitch after segmentation
img = dapi
label_img = get_img(img, model.spots, window_size=window_size, margin=math.ceil(window_size*0.1))
out = split(img, label_img, model.spots, window_size=window_size, margin=math.ceil(window_size*0.1))

# Process each tile and stitch together
cell_info={'cellid':[],'cell_center':[]}
model.spots['clustermap']=-1

for tile_split in range(out.shape[0]):
    print(f'tile: {tile_split}')
    spots_tile=out.loc[tile_split,'spots']
    dapi_tile=out.loc[tile_split,'img']

    ### instantiate model
    model_tile = ClusterMap(spots=spots_tile, dapi=dapi_tile, gene_list=gene_list, num_dims=num_dims,
                    xy_radius=xy_radius, z_radius=z_radius, fast_preprocess=False)

    if model_tile.spots.shape[0] < min_spot_per_cell:
        print(f"Less than {min_spot_per_cell} spots found in subtile. Skipping and continuing...")
        continue
    else:

        ### preprocessing
        # model_tile.preprocess(dapi_grid_interval=dapi_grid_interval, pct_filter=pct_filter)
        if sum(model_tile.spots['is_noise'] == 0) < min_spot_per_cell:
            print(f"Less than {min_spot_per_cell} non-noisy spots found in subtile. Skipping and continuing...")
            continue
        else:

            ### segmentation
            model_tile.min_spot_per_cell=min_spot_per_cell
            model_tile.segmentation(cell_num_threshold=cell_num_threshold, dapi_grid_interval=dapi_grid_interval, add_dapi=True, use_genedis=True)
            # Check if segmentation successful
            if 'clustermap' not in model_tile.spots.columns:
                continue
            else:
                # Check unique cell centers in tile
                if len(np.unique(model_tile.spots['clustermap'])) == 0:
                    print("No unique cell centers found in the cell. Skipping stitching...")
                    continue
                elif len(np.unique(model_tile.spots['clustermap'])) == 1 and np.unique(model_tile.spots['clustermap']) == [-1]:
                    print("All cell centers found were noise. Skipping stitching...")
                    continue
                else:

                    # ### stitch tiles together
                    cell_info=model.stitch(model_tile, out, tile_split)

print("Finished analyzing all subtiles")

# If tile is completely noise, throw out
if not hasattr(model, 'all_points_cellid'):
    print("No denoised cell centers found in this tile.")
    sys.exit()

# model.spots['clustermap'] = -1
# model.preprocess(dapi_grid_interval=dapi_grid_interval, pct_filter=pct_filter)
# model.min_spot_per_cell = min_spot_per_cell # cell_num_threshold proportion to tile size?
# model.segmentation(cell_num_threshold=cell_num_threshold, dapi_grid_interval=dapi_grid_interval, add_dapi=True, use_genedis=True)
model.save_segmentation(os.path.join(opath, tile, 'spots.csv'))
cell_center_df = pd.DataFrame({'cell_barcode' : model.cellid_unique.astype(int),'column' : model.cellcenter_unique[:,1],'row': model.cellcenter_unique[:,0],'z_axis':model.cellcenter_unique[:,2]})
cell_center_df.to_csv(os.path.join(opath, tile, 'cell_center.csv'))

## save figs
# preprocessing
plt.figure(figsize=(10,10))
plt.imshow(dapi.max(axis=2), cmap=plt.cm.gray)
plt.scatter(model.spots.loc[model.spots['is_noise'] == 0, 'spot_location_1'], model.spots.loc[model.spots['is_noise'] == 0, 'spot_location_2'], s=0.5, c='g', alpha=.5)
plt.scatter(model.spots.loc[model.spots['is_noise'] == -1, 'spot_location_1'], model.spots.loc[model.spots['is_noise'] == -1, 'spot_location_2'], s=0.5, c='r', alpha=.5)
plt.savefig(os.path.join(opath, tile, 'spots_pp.png'))
plt.clf()
plt.close()

## segmentation
cell_ids = model.spots['clustermap']
cells_unique = np.unique(cell_ids)
spots_repr = np.array(model.spots[['spot_location_2', 'spot_location_1']])[cell_ids>=0]
cell_ids = cell_ids[cell_ids>=0]
cmap = np.random.rand(int(max(cell_ids)+1), 3)
fig, ax = plt.subplots(figsize=(10,10))
ax.imshow(np.zeros(dapi.max(axis=2).shape), cmap='Greys_r')
ax.scatter(spots_repr[:,1],spots_repr[:,0], c=cmap[[int(x) for x in cell_ids]], s=1, alpha=.5)
ax.scatter(model.cellcenter_unique[:,1], model.cellcenter_unique[:,0], c='r', s=3)
plt.axis('off')
plt.savefig(os.path.join(opath, tile, 'cell_seg.png'))
plt.clf()
plt.close()

fig, ax = plt.subplots(figsize=(10,10))
ax.imshow(dapi.max(axis=2), cmap='Greys_r')
ax.scatter(spots_repr[:,1],spots_repr[:,0], c=cmap[[int(x) for x in cell_ids]], s=1, alpha=.5)
ax.scatter(model.cellcenter_unique[:,1], model.cellcenter_unique[:,0], c='r', s=3)
plt.axis('off')
plt.savefig(os.path.join(opath, tile, 'cell_seg_with_dapi.png'))
plt.clf()
plt.close()

stop = timeit.default_timer()
computation_time = round((stop - start) / 60, 2)

### save log (csv)
log_dict = {'number_of_spots': spots.shape[0],
           'number_of_spots_after_pp': model.spots.loc[model.spots['is_noise'] == 0, :].shape[0],
           'number_of_cells': model.cellcenter_unique.shape[0],
           'computation_time': computation_time}
log = pd.DataFrame(log_dict, index=[tile])
log.to_csv(os.path.join(opath, tile, 'log.csv'))
