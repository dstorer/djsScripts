#!/bin/bash

while getopts ":f:r:b:e:" option
do
  case $option in
    # A directory where FHD outputs are
    f) fhd_path="$OPTARG";;
    # A directory where the raw visibilities are stored
    r) raw_path="$OPTARG";;
    # Baselines to read/write
    b) baselines=$OPTARG;;
    # Function(s) to execute
    e) exec_func=$OPTARG;;
    \?) echo "Unknown option: Accepted flags are -f (obs_file_name),"
        exit 1;;
    :) echo "Missing option argument for input flag"
       exit 1;;
  esac
done

qsub -v fhd_path=${fhd_path},raw_path=${raw_path},baselines=${baselines},exec_func=${exec_func} -q hera -j oe -o ${fhd_path}/vis_read.out -l nodes=1:ppn=1 -l vmem=164G -N vis_analysis /lustre/aoc/projects/hera/dstorer/Projects/scripts/readFHD/call_FHD_IO.sh
