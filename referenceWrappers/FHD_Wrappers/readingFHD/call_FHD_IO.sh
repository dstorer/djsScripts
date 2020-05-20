#!/bin/sh

echo "JOB START TIME" `date +"%Y-%m-%d_%H:%M:%S"`
echo "Processing...."

conda activate hera

python /lustre/aoc/projects/hera/dstorer/Projects/ssinsFHDRun/beamTesting/flagged30sig/FHD_IO.py ${fhd_path} ${raw_path} ${baselines} ${exec_func}


echo "ALL DONE"
echo "JOB END TIME" `date +"%Y-%m-%d_%H:%M:%S"`
