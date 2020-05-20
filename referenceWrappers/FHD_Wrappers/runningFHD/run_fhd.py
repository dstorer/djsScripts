import numpy as np
import argparse
import os
from pyuvdata import UVData
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('obs_file', help='The path to a txt file containing the path to all uvfits files to be executed on')
parser.add_argument('outdir', help='Output directory')
args = parser.parse_args()

print('obs_file is:')
print(str(args.obs_file))
print('outdir is:')
print(str(args.outdir))

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

f = open(args.obs_file, "r")
file_names = f.read().split('\n')
print(file_names[:])

for path in file_names[:]:
	print(path)
	uv = UVData()
	uv.read(path)
	uv.reorder_blts(order='time')
	unique = np.unique(uv.time_array)
	i = 0
	while i < 21:
		print(i)
		times = unique[i:i+3]
		i = i+3
		print(times[0])
		version = str(times[0]) + '_H3C_unflagged_DipoleBeam'
		uv2 = uv.select(times=times, inplace=False)
		uv2.phase_to_time(np.median(unique))
		outname = '/lustre/aoc/projects/hera/dstorer/H3C_data/2458801/miniFiles/' + str(times[0]) + '.uvfits'
		uv2.write_uvfits(outname, spoof_nonessential=True)
		os.system("idl -e test_fhd_fromshellscript -args " + outname + " " + version + " " + args.outdir)
