#!/bin/bash


while getopts ":f:o:p:" option
do
  case $option in
    # A text file where each line is an obsid
    f) obs_file_name="$OPTARG";;
    # The output directory for the error log
    o) outdir=$OPTARG;;
    p) order=$OPTARG;;
    \?) echo "Unknown option: Accepted flags are -f (obs_file_name),"
        exit 1;;
    :) echo "Missing option argument for input flag"
       exit 1;;
  esac
done

#Manual shift to the next flag.
shift $(($OPTIND - 1))

#Throw error if no obs_id file.
if [ -z ${obs_file_name} ]; then
   echo "Need to specify a full filepath to a list of viable datapaths."
   exit 1
fi

#Throw error if no output directory
if [ -z ${outdir} ]; then
    echo "Need to specify an output directory for the error log"
fi
N_obs=$(wc -l < $obs_file_name)

qsub -v obs_file_name=${obs_file_name},N_obs=${N_obs},outdir=${outdir},order=${order}, -q hera -j oe -o ${outdir}/SSINS.out -l nodes=1:ppn=1 -l vmem=64G -N SSINS_Flagging /lustre/aoc/projects/hera/dstorer/SSINS_Scripts/call_SSINS.sh
