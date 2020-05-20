import numpy as np

import argparse
import os
import os.path
from os import path
from pyuvdata import UVData
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('obs_file', help='The path to a txt file containing the path to all uvfits files to be executed on')
parser.add_argument('outdir', help='Output directory')
parser.add_argument('version_str', help='A string to include in the name of all outputs indicating the run version')
parser.add_argument('PBS_ARRAYID', help='The index of the array job to run')
args = parser.parse_args()

print('obs_file is:')
print(str(args.obs_file))
print('outdir is:')
print(str(args.outdir))
print('Array ID is:')
print(str(args.PBS_ARRAYID))

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

ind = int(args.PBS_ARRAYID)-1
f = open(args.obs_file, "r")
file_names = f.read().split('\n')
filepath = file_names[ind]
print(filepath)

uv = UVData()
uv.read(filepath)
uv.reorder_blts(order='time')
unique = np.unique(uv.time_array)
phaseCenter = np.median(unique)
i = 0
while i < 60:
	print(i)
	times = unique[i:i+5]
	print(times[0])
	version = str(times[0]) + '_' + str(ind+1) + '_' + str(i) + '_' + args.version_str
	i = i+5
	uv2 = uv.select(times=times, inplace=False)
	uv2.phase_to_time(phaseCenter)
	outname = args.outdir + '/miniFiles/' + str(times[0]) + '.uvfits'
	dirname = args.outdir + '/fhd_' + version + '/calibration'
	if path.exists(dirname):
		print('%s has already been run - SKIPPING' % version)
		continue
	uv2.write_uvfits(outname, spoof_nonessential=True)
	os.system("idl -e run_h1c_dipolebeam_newParams -args " + outname + " " + version + " " + args.outdir)
