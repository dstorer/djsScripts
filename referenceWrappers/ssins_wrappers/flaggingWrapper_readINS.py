
from SSINS import INS, SS, MF
from SSINS import Catalog_Plot as cp
import numpy as np
import argparse
import os
from pyuvdata import UVData
from astropy.coordinates import Angle
import json

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

saved_path = '/lustre/aoc/projects/hera/iware/CHAMP/ElCap_Desktop/RFI_project/SSINS/IW_outputs/hera_SSINS_outputs'

shape_dict = {'dig1' : [1.125e8, 1.15625e8],
        'dig2' : [1.375e8, 1.40625e8],
        'dig3' : [1.625e8, 1.65625e8],
        'dig4' : [1.875e8, 1.90625e8],
        'TV4' : [1.74e8, 1.82e8],
        'TV5' : [1.82e8, 1.9e8],
        'TV6' : [1.9e8, 1.98e8],
        'TV7' : [1.98e8, 2.06e8]}

for path in file_names[:]:
	fullname = path[43:]
	name = path[43:-5]
	JD = path[47:-8]
	saved_file = saved_path + '/' + name + '.uvh5_SSINS_data.h5'
	print('Running SSINS on ' + fullname)
	#ss = SS()
	#infile = '%s%s' % (prefix, name)
	infile = path
	#ss.read(infile, flag_choice='original')
	ins = INS(saved_file)
	ang = Angle(ins.lst_array[0], unit='rad')
	ang_hour = ang.hour
	xticks = np.arange(0, len(ins.freq_array), 115)
	xticklabels = ['%.1f' % (ins.freq_array[tick]* 10 ** (-6)) for tick in xticks]
	cp.INS_plot(ins, '%s/%s_RAW' % (args.outdir, name), title='JD_%s_LST_%s' % (JD, ang_hour), vmin=0, vmax=0.03, ms_vmin=-5, ms_vmax=5, sample_sig_vmax=25,xticks=xticks,xticklabels=xticklabels,sample_sig_vmin=-25)
	ins.order=int(args.order)
	ins.metric_ms = ins.mean_subtract()
	ins.metric_array[:, :82] = np.ma.masked
	ins.metric_ms = ins.mean_subtract()

	mf = MF(ins.freq_array, 5, shape_dict=shape_dict, N_samp_thresh=25, streak=True)
	mf.apply_match_test(ins, apply_samp_thresh=True)

	cp.INS_plot(ins, '%s/%s' % (args.outdir, name), title='JD_%s_LST_%s' % (JD, ang_hour), vmin=0, vmax=0.03, ms_vmin=-5, ms_vmax=5, sample_sig_vmax=25,xticks=xticks,xticklabels=xticklabels,sample_sig_vmin=-25)


	#ss.apply_flags(flag_choice='INS', INS=ins)
	name = name + '_streakON'
	ins.write('%s/%s' % (args.outdir, name), output_type='match_events')
	ins.write('%s/%s' % (args.outdir, name), output_type='flags')
	ins.write('%s/%s' % (args.outdir, name), output_type='mask')
	#ss.write(filename_out='%s/%s_SSINSflagged.uvfits' % (args.outdir, name), file_type_out='uvfits', filename_in=infile, nsample_default=16, write_kwargs={'force_phase' : True, 'spoof_nonessential' : True})
settingsFile = open('%s/%s' % (args.outdir,'SSINS_Settings.txt'),"w")
settingsFile.write('################ SSINS Flagging Settings ################ \n')
settingsFile.write('%s is %s %s' % ('obs_file',args.obs_file,'\n'))
settingsFile.write('%s is %s %s' % ('order',args.order,'\n'))
settingsFile.write('Shape dict is: \n')
settingsFile.write(json.dumps(shape_dict))
settingsFile.write('\n')
settingsFile.write('INS is being read in from: ' + saved_path + ' \n')
settingsFile.write('\n')
settingsFile.write('##### INS Plot Settings ##### \n')
settingsFile.write('vmin=0 \n')
settingsFile.write('vmax=0.03 \n')
settingsFile.write('ms_vmin=-5 \n')
settingsFile.write('ms_vmax=5 \n')
settingsFile.write('sample_sig_vmax=25 \n')
settingsFile.write('sample_sig_vmin=-25 \n')
settingsFile.write('\n')
settingsFile.write('##### MF Settings ##### \n')
settingsFile.write('N_samp_thresh=25 \n')
settingsFile.write('streak=True \n')
settingsFile.write('apply_samp_thresh=True \n')
settingsFile.write('sig_thresh=5 \n')
