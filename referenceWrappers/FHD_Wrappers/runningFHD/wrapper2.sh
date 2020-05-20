#!/bin/sh

echo "JOB START TIME" `date +"%Y-%m-%d_%H:%M:%S"`

echo "Processing"

echo $IDL_PATH
source ~/.bashrc
echo $IDL_PATH
python /lustre/aoc/projects/hera/dstorer/Projects/IDLscripts/run_fhd.py ${obs_file_name} ${outdir}

echo "JOB END TIME" `date +"%Y-%m-%d_%H:%M:%S"`
