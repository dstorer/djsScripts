#!/bin/sh

echo "JOB START TIME" `date +"%Y-%m-%d_%H:%M:%S"`

echo "There are $N_obs observations in this file"

echo "Processing"


conda activate hera

python /lustre/aoc/projects/hera/dstorer/SSINS_Scripts/multiSigThreshWrapper.py ${obs_file_name} ${outdir} ${order}


echo "JOB END TIME" `date +"%Y-%m-%d_%H:%M:%S"`
