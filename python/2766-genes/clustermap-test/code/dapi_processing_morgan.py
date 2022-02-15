"""
Tile-by-tile pre-processing of DAPI stains for Peng Intestine Samples, for input into ClusterMap
Performs image pre-processing on each DAPI stain to homogenize DAPI intensities and filter out tissue noise
Performs image pre-processing of corresponding flamingo stain and creates a mask to increase resolution of DAPI boundaries
Applies flamingo mask to DAPI stain and writes to file

Inputs to the script are tile number and output folder
Outputs of the script are processed 3D DAPI tiles
"""

# Load libraries
import argparse
import sys
import tifffile as tiff
import numpy as np
import os
import time
from scipy import ndimage
from skimage.morphology import ball, opening, closing, dilation
from skimage.exposure import rescale_intensity
from skimage.filters import threshold_otsu

# ============================== DEFINE HELPER FNS FOR FILE I/O ==============================

# Load images (DAPI and Flamingo stains)
def loadStain(tile_num, stain):
    """
    Load either DAPI stain or flamingo stain of a specified tile number
    Returns image stored in NumPy arrays
    """
    data_dir = '/stanley/WangLab/morgan/PengIntestines/processed_data/output/max'
    path = os.path.join(data_dir, f'{stain}/tile_{tile_num}.tif')
    img = tiff.imread(path)

    return img

# Save masked DAPI image
def saveDapi(tile, outpath):
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    tiff.imwrite(os.path.join(outpath, f'tile_{tile.num}.tif'), tile.dapi_masked)

# ============================== DEFINE TILE CLASS ==============================

class Tile():
    """
    Variables:
        *.num: number of tile or position
		*.shape: z, col, row dimensions of image
		*.rotated: boolean indicating whether or not stains have been rotated
        *.dapi: 3D DAPI stain
        *.dapi_processed: pre-processed 3D DAPI image
        *.dapi_masked: dapi_processed masked by flamingo_mask
        *.flamingo: 3D flamingo stain
        *.flamingo_mask: pre-processed and binarized flamingo stain

    Functions:
		*.rotateImg
        *.createFlamingoMask
        *.processDapi
        *.maskDapi
    """

    def __init__(self, tile_num):
        # Constructor 1: Provide only tile number
        self.num = tile_num
        self.dapi = loadStain(self.num, 'dapi')
        self.flamingo = loadStain(self.num, 'flamingo')
        self.rotated = False

        # Check if dimensions of dapi and flamingo match
        if (self.dapi.shape == self.flamingo.shape):
            self.shape = self.dapi.shape
        else:
            print(f"Error: mismatched dimensions. DAPI dimensions are {self.dapi.shape} and flamingo dimensions are {self.flamingo.shape}")
            sys.exit()

    #def rotateImg(self, left_rotation):
    #    for z in range(self.shape[0]):
    #        self.dapi[z,:,:] = ndimage.rotate(self.dapi[z,:,:], left_rotation, reshape=False)
    #        self.flamingo[z,:,:] = ndimage.rotate(self.flamingo[z,:,:], left_rotation, reshape=False)
    #    self.rotated = True

    def createFlamingoMask(self, strel = ball(3), ub = 90, lb = 1):
        start = time.time()
        # Rescale intensity between upper and lower bounds of intensity percentiles
        p_lower, p_upper = np.percentile(self.flamingo, (lb, ub))
        rescaled = rescale_intensity(self.flamingo, in_range=(p_lower, p_upper))
        # Morphological opening
        opened = opening(rescaled, strel)
        # Otsu thresholding
        self.flamingo_mask = opened > threshold_otsu(opened)
        print(f"Created flamingo mask [{round(time.time() - start,3)}s]")

    def processDapi(self, strel=ball(3), filter_bg_threshold=95):
        """
        Pre-processing on DAPI:
        Sequential operations are performed:
            1) Intensity Filtering
            2) Morphological closing for homogenizing intensity of DAPI stain
            3) Dilation
        """
        start = time.time()
        # Intensity filter original dapi stain
        filter_value = np.percentile(self.dapi[self.dapi != 0], filter_bg_threshold)
        dapi_filtd = np.where(self.dapi > filter_value, self.dapi, 0)
        # Morphological closing
        dapi_filtd_closed = closing(dapi_filtd, strel)
        # Dilation
        self.dapi_processed = dilation(dapi_filtd_closed, strel)
        print(f"Finished DAPI processing [{round(time.time() - start, 3)}s]")

    def maskDapi(self):
        """
        Apply flamingo mask to processed DAPI
        """
        self.dapi_masked = np.where(self.flamingo_mask == True, 0, self.dapi_processed)

# ============================== MAIN ==============================

def main(args):
    startTime = time.time()

    # Parse arguments
    tile_num = args.tile_num
    outpath = args.outpath

	# Initialize tile object
    tile = Tile(tile_num)
    print("Tile object created")

	# Perform processing
    tile.createFlamingoMask()
    tile.processDapi()
    tile.maskDapi()

	# Save masked DAPI
    saveDapi(tile, outpath)
    print(f"Tile {tile_num} masked dapi saved to file [Total time: {round(time.time() - startTime, 3)}s]")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Tile-by-tile preprocessing and masking of DAPI stain using flamingo stain")
    parser.add_argument('tile_num', type=str, help='tile number to be processed')
    parser.add_argument('outpath', type=str, help='path to dapi_masked output directory')
    args = parser.parse_args()
    main(args)
