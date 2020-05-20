from __future__ import print_function

from pyuvdata import UVData
from pyuvdata import UVCal
import os
import subprocess
import os.path
from os import path
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
from matplotlib import pyplot as plt
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, Angle
from astropy.time import Time
import matplotlib.animation as anm
from scipy.interpolate import interp1d
import helperFunctions
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('fhd_path', help='The path to a directory containing the desired FHD outputs')
parser.add_argument('raw_path')
parser.add_argument('baselines', help='A list of tuples of baselines to read and store visibilities for')
parser.add_argument('exec_func')
args = parser.parse_args()

fhd_path = args.fhd_path
baselines = args.baselines
rawVisPath = args.raw_path
exec_func = args.exec_func
#freqrange=[120,185] #Frequency range (MHz) to read in and plot
#visrange=[1,-1] #Range (indices) of observations array to read in and plot
#Option to write out a calfits file:
write_calfits = False

#### Setup paths and arrays ####

thedir = fhd_path
vis = [name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name))]
vis.sort()
vis = vis[1:-1]
obspath = fhd_path + '/' + vis[4]
data_array = {}

#Get telescope location:
dat = UVData()
version = vis[4]
print(version)
obsid = version.split('_')[1]
prefix = obspath + '/' + 'metadata' + '/' + obsid + '_'
files1 = [prefix + f for f in ['params.sav','settings.txt', 'layout.sav']]
prefix = obspath + '/' + 'vis_data' + '/' + obsid + '_'
files2 = [prefix + f for f in ['flags.sav','vis_XX.sav','vis_YY.sav','vis_model_XX.sav','vis_model_YY.sav']]
files = np.append(files1, files2)

dat.read_fhd(files)
loc = EarthLocation.from_geocentric(*dat.telescope_location, unit='m')

print('rawVisPath is %s' % rawVisPath)


if exec_func == 'all' or 'raw' in exec_func:
    print('#################################')
    print('Reading in raw visibilities......')
    print('#################################')

    vis_array = helperFunctions.readRawVisibilities(dat,path1=rawVisPath, baselines=baselines)
    with open('%s/raw_vis_array.pickle' % fhd_path, 'wb') as f:
        pickle.dump(vis_array, f)
    del vis_array
    
if exec_func == 'all' or 'dirty' in exec_func:
    print('#######################################')
    print('Reading in FHD dirty visibilities......')
    print('#######################################')

    vis_array = helperFunctions.readFHDVisibilities(dat,path1=fhd_path, baselines=baselines, use_model=False)
    with open('%s/fhd_dirty_vis.pickle' % fhd_path, 'wb') as f:
        pickle.dump(vis_array, f)
    del vis_array
    
if exec_func == 'all' or 'model' in exec_func:
    print('#######################################')
    print('Reading in FHD model visibilities......')
    print('#######################################')

    vis_array = helperFunctions.readFHDVisibilities(dat,path1=fhd_path, baselines=baselines, use_model=True)
    with open('%s/fhd_model_vis.pickle' % fhd_path, 'wb') as f:
        pickle.dump(vis_array, f)
    del vis_array

if exec_func == 'all' or 'cal' in exec_func:
    print('#######################################')
    print('Reading in FHD model visibilities......')
    print('#######################################')

    vis_array = helperFunctions.readCalSolutions(dat,path1=fhd_path, metric='amp')
    with open('%s/fhd_gains_amp.pickle' % fhd_path, 'wb') as f:
        pickle.dump(vis_array, f)
    del vis_array

    vis_array = helperFunctions.readCalSolutions(dat,path1=fhd_path, metric='phase')
    with open('%s/fhd_gains_phase.pickle' % fhd_path, 'wb') as f:
        pickle.dump(vis_array, f)
    del vis_array

print('ALL DONE!')

