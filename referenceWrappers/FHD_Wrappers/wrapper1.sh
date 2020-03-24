#!/bin/bash

while getopts ":f:o:p:" option
do
  case $option in
    # A text file where each line is an obsid
    f) obs_file_name="$OPTARG";;
    # The output directory for the error log
    o) outdir=$OPTARG;;
    \?) echo "Unknown option: Accepted flags are -f (obs_file_name),"
        exit 1;;
    :) echo "Missing option argument for input flag"
       exit 1;;
  esac
done
#t=${obs_file_name:54:-4}
t=${obs_file_name:87:-4}
echo ${t}
qsub -v obs_file_name=${obs_file_name},outdir=${outdir} -j oe -o ${outdir}/FHD_${t}.out -l nodes=1:ppn=1 -l vmem=128G -N FHD_cal_${t} /lustre/aoc/projects/hera/dstorer/Projects/IDLscripts/wrapper2.sh 
