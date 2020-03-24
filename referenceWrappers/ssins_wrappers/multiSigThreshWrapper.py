
from SSINS import INS, SS, MF
from SSINS import Catalog_Plot as cp
import numpy as np
import argparse
import os
from pyuvdata import UVData

parser = argparse.ArgumentParser()
parser.add_argument('obs_file', help='The path to a txt file containing a list of obsids')
parser.add_argument('outdir', help='The output directory')
parser.add_argument('order', help='The order of the INS polynomial fit')
args = parser.parse_args()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

#indir = args.infile[:args.infile.rfind('/')]
#if indir == args.outdir:
#    raise ValueError("indir and outdir are the same")

f = open(args.obs_file, "r")
file_names = f.read().split('\n')

shape_dict = {'TV4': [174e6, 182e6],
                  'TV5': [182e6, 190e6],
                  'TV6': [190e6, 198e6],
                  'TV7': [198e6, 206e6],
                  'TV8': [206e6, 214e6],
                  'TV9': [214e6, 222e6],
                  'TV10': [222e6, 230e6],
                  'TV11': [230e6, 238e6],
                  'AM1': [45e6, 57e6],
                  'AM2': [63e6, 70e6],
                  'AM3': [73e6, 80e6],
                  'MB1': [128e6, 133e6],
                  'MB2': [140e6, 147e6],
                  'MB3': [153e6, 160e6],
                  'MB4': [165e6, 170e6]}

sig_thresh = {
	'TV4' : 5,
	'TV5' : 5,
	'TV6' : 5,
	'TV7' : 5,
	'TV8' : 5,
	'TV9' : 5,
	'TV10' : 5,
	'TV11' : 5,
        'AM1': 5,
        'AM2': 5,
        'AM3': 5,
        'MB1': 5,
        'MB2': 5,
        'MB3': 5,
        'MB4': 5,
	'narrow' : 5,
	'streak' : 20
}

dab_width = 1.536e6
dab_freqs = np.arange(214e6, 230e6 + dab_width, dab_width)
dab_dict = {'DAB%i' % ind: [dab_freqs[ind], dab_freqs[ind + 1]]  for ind in range(len(dab_freqs) - 1)}

for shape in dab_dict:
	shape_dict[shape] = dab_dict[shape]
	sig_thresh[shape] = 5

for path in file_names[:]:
	fullname = path[51:]
	name = path[51:-5]
	print('Running SSINS on ' + fullname)
	ss = SS()
	#infile = '%s%s' % (prefix, name)
	infile = path
	ss.read(infile, flag_choice='original')
	ins = INS(ss)
	#cp.INS_plot(ins, '%s/%s_RAW' % (args.outdir, name), vmin=0, vmax=0.03, ms_vmin=-5, ms_vmax=50)
	#ins.order=int(args.order)
	ins.metric_ms = ins.mean_subtract()
	#ins.write('%s/%s' % (args.outdir, name))
	where_FM = np.where(np.logical_and(ins.freq_array > 87.5e6, ins.freq_array < 108e6))
	ins.metric_array[:, where_FM] = np.ma.masked
	ins.metric_ms = ins.mean_subtract()
	#ins.metric_array[:, :82] = np.ma.masked

	mf = MF(ins.freq_array, sig_thresh, shape_dict=shape_dict, N_samp_thresh=10, streak=True)
	mf.apply_match_test(ins, apply_samp_thresh=True)

	#cp.INS_plot(ins, '%s/%s' % (args.outdir, name), vmin=0, vmax=0.03, ms_vmin=-5, ms_vmax=5)


	ss.apply_flags(flag_choice='INS', INS=ins)
	name = name + '_streak20sig'
	#ins.write('%s/%s' % (args.outdir, name), output_type='match_events')
	#ins.write('%s/%s' % (args.outdir, name), output_type='flags')
	#ins.write('%s/%s' % (args.outdir, name), output_type='mask')
	ss.write(filename_out='%s/%s_SSINSflagged_streak20.uvh5' % (args.outdir, name), file_type_out='uvh5', filename_in=infile, nsample_default=16)
