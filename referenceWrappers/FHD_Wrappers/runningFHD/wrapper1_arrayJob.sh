#!/bin/bash

while getopts ":f:o:v:n:" option
do
  case $option in
    # A text file where each line is an obsid
    f) obs_file_name="$OPTARG";;
    # The output directory for the error log
    o) outdir="$OPTARG";;
    # The run version
    v) version_str=$OPTARG;;
    # A number or prefix to include in the name of the output log
    n) num_prefix=$OPTARG;;
    \?) echo "Unknown option: Accepted flags are -f (obs_file_name),"
        exit 1;;
    :) echo "Missing option argument for input flag"
       exit 1;;
  esac
done
#t=${obs_file_name:54:-4}
t=${obs_file_name:72:-4}
h=${obs_file_name:87:-4}
echo ${t}
qsub -v obs_file_name=${obs_file_name},outdir=${outdir},version_str=${version_str} -j oe -o ${outdir}/FHD_${num_prefix}.out -l nodes=1:ppn=1 -l mem=64G -l vmem=64G -N FHD_cal_${num_prefix} -q hera /lustre/aoc/projects/hera/dstorer/Projects/IDLscripts/standardFhdRun/wrapper2_arrayJob.sh 
