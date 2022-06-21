#!/bin/bash -l
ppath="/stanley/WangLab/jiahao/Github/mAD-analysis/python/2766-genes/clustermap-test"
#$ -o /stanley/WangLab/jiahao/Github/mAD-analysis/python/2766-genes/clustermap-test/log_o/qsub_log_o.$JOB_ID.$TASK_ID
#$ -e /stanley/WangLab/jiahao/Github/mAD-analysis/python/2766-genes/clustermap-test/log_e/qsub_log_e.$JOB_ID.$TASK_ID

source "/broad/software/scripts/useuse"
reuse Anaconda3
source activate /stanley/WangLab/envs/clustermap
now=$(date +"%T")
echo "Current time : $now"
position=$(sed "${SGE_TASK_ID}q;d" /stanley/WangLab/jiahao/Github/mAD-analysis/python/2766-genes/clustermap-test/code/position.txt)
position=$(echo $position | tr -d '\r')
python $ppath/code/mad_clustermap_test.py $position

echo "Finished"
now=$(date +"%T")
echo "Current time : $now"
