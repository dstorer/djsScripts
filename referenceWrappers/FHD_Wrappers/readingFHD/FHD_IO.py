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
#rawVisPath = "/lustre/aoc/projects/hera/H1C_IDR2/IDR2_2/2458098"
#JDs to execute on (assumes FHD results are organized into sub-directories by these names within path1)
days = ['2458098']
#ant_nums = [51,52] #Antenna numbers to read in
#sampleAnt = 51
#baselines = [(51,51),(70,70),(52,52),(53,53)]
#plotName = 'H3C_SSINSFlagged_StreakOn' #Base name for plots
#outpath = "/lustre/aoc/projects/hera/dstorer/H3C_data/2458787/FHD_Output/streakOn/plots/" #Output directory for plots
pol = 'Jxx' #Polarization to read in
freqrange=[120,185] #Frequency range (MHz) to read in and plot
visrange=[1,-1] #Range (indices) of observations array to read in and plot
colormap = 'plasma_r'
#Source to read calibration from (Either 'sav' for FHD .sav files, or 'fits' for UVCal cal.fits files:
#Note that for cal.fits files I assume they are stored in the calibration subdirectory with the same naming convention
source = 'sav'
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
#print(obsid)
prefix = obspath + '/' + 'metadata' + '/' + obsid + '_'
files1 = [prefix + f for f in ['params.sav','settings.txt', 'layout.sav']]
prefix = obspath + '/' + 'vis_data' + '/' + obsid + '_'
#print(prefix)
files2 = [prefix + f for f in ['flags.sav','vis_XX.sav','vis_YY.sav','vis_model_XX.sav','vis_model_YY.sav']]
files = np.append(files1, files2)
#print(files)

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

