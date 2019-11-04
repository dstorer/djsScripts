from pyuvdata import UVData
from pyuvdata import UVCal
import os
import os.path
from os import path
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, Angle
from astropy.time import Time

#Path to FHD run:
path1 = "/lustre/aoc/projects/hera/dstorer/Projects/ssinsFHDRun/results/"
#days = ['2458098','2458099','2458101','2458102']
days = ['2458098']
#days = ['2458098','2458099']
ant_nums = [51,70,52,53]
outpath = "/lustre/aoc/projects/hera/dstorer/Projects/ssinsFHDRun/results/plots/"
plotName = 'masked_calibrationSolutionWaterfall_ALL.png'
thedir = path1 + days[0]
vis = [ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]
#vis = filter(os.path.isdir, os.listdir(path1 + days[0]))
print(len(vis))
print(vis[7])
vis.sort()
obspath = path1 + days[0] + '/' + vis[1]
pol = 'Jxx'
#
#Source to read calibration from (Either 'sav' for FHD .sav files, or 'fits' for UVCal cal.fits files:
#Note that for cal.fits files I assume they are stored in the calibration subdirectory with the same naming convention
source = 'sav' 
#Option to write out a calfits file:
write_calfits = False
freqrange=[120,180]
colormap = 'plasma_r'
mask = True
threshold = 0.0005

#################

#Get telescope location:
dat = UVData()
version = vis[1]
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

nants = len(ant_nums)
for k in range(len(days)):
    thedir = path1 + days[k]
    vis = [ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]
    vis.sort()
    fig, axs = plt.subplots(1, nants, figsize=(15,6))
    fig.subplots_adjust(hspace=0.5,wspace=0.5)
    fig.subplots_adjust(right=0.8)
    for j in range(nants):
	nocal = 0
        i=0
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
                freq_array = cal.freq_array[0,:]
                freq_array = np.transpose(freq_array)
                freq_array = np.divide(freq_array,1000000)
                minfreqind = np.argmax(freq_array>freqrange[0])
                maxfreqind = np.argmin(freq_array<freqrange[1])
		print('Resetting cal_array')
                cal_array = np.zeros((nfreqs,nobs))
		mask_array = np.zeros((nfreqs,nobs))
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
		    mask_array[:,i] = np.ones(np.shape(gain))[:,0]
		    i = i + 1
		    continue
            cal_array[:,i] = gain[:,0]
	    i = i + 1
        cal_array = np.transpose(cal_array)
	mask_array = np.transpose(mask_array)
        cal_array = cal_array[:,minfreqind:maxfreqind]
	mask_array = mask_array[:,minfreqind:maxfreqind]
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
	## Getting RA of zenith ##
	zenith_start = SkyCoord(Angle(0, unit='deg'),Angle(90,unit='deg'),frame='altaz',obstime=obstime_start,location=loc)
	zenith_start = zenith_start.transform_to('icrs')
	zenith_end = SkyCoord(Angle(0, unit='deg'),Angle(90,unit='deg'),frame='altaz',obstime=obstime_end,location=loc)
	zenith_end = zenith_end.transform_to('icrs')
	print(zenith_start.ra.degree)
	print(zenith_end.ra.to('rad'))
        ## Plotting ##
	mx = np.ma.masked_array(cal_array, mask=mask_array)
	currmap = mpl.cm.get_cmap()
	currmap.set_bad('black',1.)
        im = axs[j].imshow(
            mx, origin='upper', interpolation='none',
            cmap = colormap,
            extent=[
            freq_array[minfreqind], freq_array[maxfreqind],
            endTime, startTime
            ],
            aspect='auto', vmin=0, vmax=0.2
        )
	JD = int(obstime_start.jd)
        axs[j].set_title('Antenna #' + str(ant_nums[j]) + ', JD: ' + str(JD))
        axs[j].set_xlabel('Frequency (MHz)')
	if j == 0:
            axs[j].set_ylabel('Time (LST)')
	elif j == (nants-1):
	    axs[j].set_ylabel('RA of zenith (deg)',rotation=270,labelpad=20)
	    print(float(zenith_end.ra.degree))
	    axs[j].set_yticklabels(np.around(np.linspace(zenith_end.ra.degree,zenith_start.ra.degree,9),1))
	    axs[j].yaxis.tick_right()
	    axs[j].yaxis.set_label_position("right") 
	else:
	    axs[j].get_yaxis().set_ticklabels([])
	print('###############################################################')
	print('Out of ' + str(len(vis)) + ' observations, ' + str(nocal) + ' did not contain calibration solutions')
	print('###############################################################')
    cbar_ax = fig.add_axes([0.87,0.1,0.02,0.8])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Gain Amplitude', rotation=270, labelpad=20) 
    fig.savefig(outpath + str(days[k]) + '_' + plotName)
    print('Saving figure to: ' + outpath + str(days[k]) + '_' + plotName)
print('###############################################################')
print('##################### FINISHED RUNNING ########################')
print('###############################################################')

