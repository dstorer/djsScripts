from __future__ import print_function

from pyuvdata import UVData
from pyuvdata import UVCal
import os
import subprocess
import os.path
from os import path
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, Angle
from astropy.time import Time
import matplotlib.animation as anm
from scipy.interpolate import interp1d

#Path to FHD run:
path1 = "/lustre/aoc/projects/hera/dstorer/Projects/ssinsFHDRun/results/"
#days = ['2458098','2458099','2458101','2458102']
days = ['2458098','2458099']
ant_nums = [51,70,52,53]
outpath = "/lustre/aoc/projects/hera/dstorer/Projects/ssinsFHDRun/results/plots/"
thedir = path1 + days[0]
vis = [ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]
#vis = os.listdir(path1 + days[0])
vis.sort()
obspath = path1 + days[0] + '/' + vis[10]
pol = 'Jxx'
frame_baseName = 'calMovieFrames'
#frame_baseName = 'testNonzeroFailure'
#
#Source to read calibration from (Either 'sav' for FHD .sav files, or 'fits' for UVCal cal.fits files:
#Note that for cal.fits files I assume they are stored in the calibration subdirectory with the same naming convention
source = 'sav' 
#Option to write out a calfits file:
write_calfits = False
freqrange=[150,185]
colormap = 'plasma_r'
mask = True
threshold = 0.0005

#################

#Get telescope location:
dat = UVData()
version = vis[10]
obsid = version[4:-20]
prefix = obspath + '/' + 'metadata' + '/' + obsid + '_'
files1 = [prefix + f for f in ['params.sav','settings.txt']]
prefix = obspath + '/' + 'vis_data' + '/' + obsid + '_'
files2 = [prefix + f for f in ['flags.sav','vis_XX.sav','vis_YY.sav','vis_model_XX.sav','vis_model_YY.sav']]
files = np.append(files1, files2)

dat.read_fhd(files)
loc = EarthLocation.from_geocentric(*dat.telescope_location, unit='m')

###################

cal = UVCal()
data_array = []
nants = len(ant_nums)
for k in range(len(days)):
    data_array = []
    frame_baseName = str(days[k]) + '_' + frame_baseName
    thedir = path1 + days[k]
    vis = [ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]
    #vis = os.listdir(path1 + days[k])
    vis.sort()
    fig, axs = plt.subplots(figsize=(15,15))
    #fig.subplots_adjust(hspace=0.5,wspace=0.5)
    #fig.subplots_adjust(right=0.8)
    for j in range(nants):
	nocal = 0
        i=0
	times = []
        for v in vis:
	    data = True
	    obs = v[4:-20]
            print('Reading data for antenna ' + str(ant_nums[j]) + ', Obs ' + str(obs))
            path = path1 + days[k] + '/' + v
	    if os.path.exists(path + '/calibration/') is False:
		print('!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!')
		print(' There is no calibration data for this observation. The gains for this observation will be entered as 0.')
		nocal = nocal + 1
		data = False
		if not i == 0:
		    times.append(times[-1]+0.0139)
		else:
		    times.append(0)
            if source == 'sav':
		if data == True:
		    obsfile = path + '/metadata/' + obs + '_obs.sav'
		    calfile = path + '/calibration/' + obs + '_cal.sav'
                    settingsfile = path + '/metadata/' + obs + '_settings.txt'
                    cal.read_fhd_cal(calfile,obsfile,settings_file=settingsfile)
            elif source == 'fits':
                calfits = path + str(nums[i]) + '/calibration/' + obs + '_cal.fits'
                cal.read_calfits(calfits)
            else:
                print('####################### ERROR #####################')
		print('The source must be set to either sav or fits')

            if write_calfits == True and data == True:
                outpath = path + str(nums[i]) + '/calibration/' + obs + '_cal.fits'
                cal.write_calfits(outpath, clobber=True)
            if i == 0 and data==True:
                nobs = len(vis)
                nfreqs = cal.Nfreqs
		print('nfreqs is: ' + str(nfreqs))
                print(np.shape(cal.freq_array))
		freq_array = cal.freq_array[0,:]
                freq_array = np.transpose(freq_array)
                freq_array = np.divide(freq_array,1000000)
                minfreqind = np.argmax(freq_array>freqrange[0])
                maxfreqind = np.argmin(freq_array<freqrange[1])
                cal_array = np.empty((nfreqs,nobs))
                print(np.shape(cal_array))
		if j == 0:
                    ## Get the LST start and end times for this obsid ##
                    time_array = cal.time_array
		    obstime_start = Time(time_array[0],format='jd',location=loc)
                    startTime = obstime_start.sidereal_time('mean').hour
                    JD = int(obstime_start.jd)
		    JDtest = int(time_array[0])
		    print('JD test value is: ' + str(JDtest))
		    print('The Julian Date of this obs is: ' + str(JD))
	    if v == vis[-1] and j == 0:
		time_array = cal.time_array
		obstime_end = Time(time_array[-1],format='jd',location=loc)
		endTime = obstime_end.sidereal_time('mean').hour
	    if data is True:
                gain = np.abs(cal.get_gains(ant_nums[j], pol))
	    if data is False:
		if i == 0:
		    continue
		else:
		    gain = np.zeros(np.shape(gain))
		    print('Zeroing gains for failed calibration')
		    cal_array[:,i] = gain[:,0]
		    i = i + 1
		    continue
            cal_array[:,i] = gain[:,0]
	    i = i + 1
	    tArr = cal.time_array
	    t = Time(tArr[0],format='jd',location=loc).sidereal_time('mean').hour
	    times.append(t)
        cal_array = np.transpose(cal_array)
        cal_array = cal_array[:,minfreqind:maxfreqind]
        ## Getting max and min values ##
        calsort = np.sort(cal_array)
        maxind = int(len(calsort[0])/15)
        maxuse = calsort[:,-maxind]
        maxuse = np.max(maxuse)
	if maxuse >=5:
		maxind = int(len(calsort[0])/12)
		maxuse = calsort[:,-maxind]
		maxuse = np.max(maxuse)
		print('!!!!!!!!!!!!!!!!!!! Had to adjust max used value !!!!!!!!!!!!!!!!!!!')
        maxval = np.max(calsort[:,-1])
        print('Maximum gain value is: ' + str(maxval))
        print('Maximum gain value used is: ' + str(maxuse))
	freq_array = freq_array[minfreqind:maxfreqind]
	res = {
	        "ant_num": ant_nums[j],
	        "obsid": obs,
	        "cal_array": cal_array,
	        "time_array": times,
		"freq_array": freq_array
	    }
        data_array.append(res)
    ## Plotting ##
    fig = plt.figure()
    for t in range(len(vis)):
	for m in range(nants):
	    lbl = 'Antenna ' + str(ant_nums[m])
	    plt.scatter(freq_array,data_array[m]['cal_array'][t,:],s=2,label=lbl)
	fname = frame_baseName + '_' + str(t) + '.png'
	plt.title('Calibration Solutions for ' + str(JD))
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Gain Amplitude')
	plt.ylim(0,0.3)
	plt.legend(loc='upper right')
	ax = plt.gca()
	currTime = data_array[m]['time_array'][t]
	timeStr = 'LST: ' + str(np.around(currTime,2))
	plt.text(0.03,0.9,timeStr,transform=ax.transAxes,fontsize=20)
	plt.savefig(fname)
	plt.cla()
print('###############################################################')
print('##################### FINISHED RUNNING ########################')
print('###############################################################')

