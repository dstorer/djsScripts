#!/bin/sh

echo "JOB START TIME" `date +"%Y-%m-%d_%H:%M:%S"`
echo "Processing...."

conda activate hera

python /lustre/aoc/projects/hera/dstorer/Projects/scripts/readFHD/FHD_IO.py ${fhd_path} ${raw_path} ${baselines} ${exec_func}


echo "ALL DONE"
echo "JOB END TIME" `date +"%Y-%m-%d_%H:%M:%S"`
