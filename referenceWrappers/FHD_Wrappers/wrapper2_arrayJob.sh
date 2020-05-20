#!/bin/sh
#PBS -t 1%12
#PBS -l hostlist="herapost[001-003]"

echo "JOB START TIME" `date +"%Y-%m-%d_%H:%M:%S"`

echo "Processing"

echo $IDL_PATH
source ~/.bashrc
echo $IDL_PATH
echo ${PBS_ARRAYID}

python /lustre/aoc/projects/hera/dstorer/Projects/IDLscripts/standardFhdRun/run_fhd_h1c_arrayJob.py ${obs_file_name} ${outdir} ${version_str} ${PBS_ARRAYID}

echo "JOB END TIME" `date +"%Y-%m-%d_%H:%M:%S"`
