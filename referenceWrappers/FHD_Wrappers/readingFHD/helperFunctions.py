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

def get_zenithRA(time_array,telescope_location):
    zens = []
    for t in time_array:
        if t is None:
            zens.append(zens[-1] + (zens[-1]-zens[-2]))
        else:
            zen = SkyCoord(Angle(0, unit='deg'),Angle(90,unit='deg'),frame='altaz',obstime=t,location=telescope_location)
            zen = zen.transform_to('icrs')
            zen = zen.ra.degree
            zens.append(zen)
    return zens

def get_LSTs(time_array):
    lst_arr = []
    for i in range(len(time_array)):
        t = time_array[i]
        if t is not None:
            lst = t.sidereal_time('mean').hour
            lst_arr.append(lst)
        else:
            if i == 0:
                lst = 0
            else:
                lst = lst_arr[-1]+0.0139
            lst_arr.append(lst)
    return lst_arr

def readFHDVisibilities(uv, path1, baselines='all', days=None, polarization='XX', use_model=False, visrange='all'):
    vis_array = {'amp': {}, 'phase': {}}
    if baselines=='all':
        baselines = uv.get_antpairs()
    if days == None:
        days = [str(uv.time_array[0]).split('.')[0]]
        print(days)
    for k in range(len(days)):
        vis = [ name for name in os.listdir(path1) if os.path.isdir(os.path.join(path1, name)) and 'fhd_' in name ]
        vis.sort()
        if visrange != 'all':
            vis = vis[visrange[0]:visrange[1]]
        obs1 = vis[3].split('_')[1]
        path = '%s/%s' % (path1,vis[3])
        uv = UVData()
        fhd_files = []
        fhd_files.append(path + '/metadata/' + obs1 + '_params.sav')
        fhd_files.append(path + '/metadata/' + obs1 + '_settings.txt')
        fhd_files.append(path + '/metadata/' + obs1 + '_layout.sav')
        vis_files = ['flags.sav','vis_XX.sav','vis_YY.sav','vis_model_XX.sav','vis_model_YY.sav']
        for f in vis_files:
            fhd_files.append(path + '/vis_data/' + obs1 + '_' + f)
        uv.read_fhd(fhd_files, use_model=use_model)
        nfreqs = len(uv.freq_array[0])
        ntimes = len(np.unique(uv.time_array))
        freq_array = uv.freq_array[0]*10**(-6)
        day_array = {}
        day_array_phase = {}
        for b in baselines:
            day_array[b] = np.zeros((len(vis)*ntimes,nfreqs))
            day_array_phase[b] = np.zeros((len(vis)*ntimes,nfreqs))
        for i in range(len(vis)):
            v = vis[i]
            obs = v.split('_')[1]
            print(obs)
            path = '%s/%s' % (path1,v)
            if os.path.exists(path + '/vis_data/') is False:
                if i == 0:
                    continue
                mask_array[i*ntimes:(i+1)*ntimes,:] = np.ones(np.shape(vis_data))
                for b in baselines:
                    day_array[b][i*ntimes:(i+1)*ntimes,:] = vis_data
                    day_array_phase[b][i*ntimes:(i+1)*ntimes,:] = vis_data_phase
            else:
                uv = UVData()
                fhd_files = []
                fhd_files.append(path + '/metadata/' + obs + '_params.sav')
                fhd_files.append(path + '/metadata/' + obs + '_settings.txt')
                fhd_files.append(path + '/metadata/' + obs + '_layout.sav')
                vis_files = ['flags.sav','vis_XX.sav','vis_YY.sav','vis_model_XX.sav','vis_model_YY.sav']
                for f in vis_files:
                    fhd_files.append(path + '/vis_data/' + obs + '_' + f)
                uv.read(fhd_files, use_model=use_model, bls=baselines, file_type='fhd')
                ntimes = len(np.unique(uv.time_array))
                if i == 0:
                    mask_array = np.zeros((ntimes*len(vis),nfreqs))
                for b in baselines:
                    vis_data = np.abs(uv.get_data(b[0], b[1], polarization))
                    vis_data_phase = np.angle(uv.get_data(b[0], b[1], polarization))
                    day_array[b][i*ntimes:(i+1)*ntimes,:] = vis_data
                    day_array_phase[b][i*ntimes:(i+1)*ntimes,:] = vis_data_phase
        for b in baselines:
            day_array[b] = np.ma.masked_array(day_array[b],mask=mask_array)
            day_array_phase[b] = np.ma.masked_array(day_array_phase[b],mask=mask_array)
        vis_array['amp'][days[k]] = day_array
        vis_array['phase'][days[k]] = day_array_phase
    return vis_array

def readRawVisibilities(uv, path1, baselines='all', days=None, polarization='XX', jdrange=None):
    vis_array = {}
    if baselines == 'all':
        baselines = uv.get_antpairs()
    if days == None:
        days = [str(uv.time_array[0]).split('.')[0]]
        print(days)
    for k in range(len(days)):
        day_array = {}
        for b in baselines:
            day_array[b] = []
        vis = [ name for name in os.listdir(path1) if name.endswith(".HH.uvh5") ]
        vis.sort()
        if jdrange != None:
            jds = [float(name.split('.')[2]) for name in vis]
            minind = np.argmin(np.abs(np.subtract(jds,jdrange[0])))
            maxind = np.argmin(np.abs(np.subtract(jds,jdrange[1])))
            vis = vis[minind:maxind]
        for i in range(len(vis)):
            v = vis[i]
            print(v)
            path = '%s/%s' % (path1,v)
            uv = UVData()
            if baselines == 'all':
                uv.read(path)
            else:
                uv.read(path, bls=baselines)
            for b in baselines:
                vis_data = uv.get_data(b[0], b[1], polarization)
                for x in vis_data[:]:
                    day_array[b].append(x)
        vis_array[days[k]] = day_array
    return vis_array

def readCalSolutions(uv,path1,data_array={},days=None,ant_nums='all',pol='xx',visrange='all',source='sav',
              write_calfits=False, metric='amp'):
    cal = UVCal()
    if ant_nums == 'all':
        ant_nums = uv.antenna_numbers
    ant_nums.sort()
    nants = len(ant_nums)
    if days == None:
        days = [str(uv.time_array[0]).split('.')[0]]
        print(days)
    for k in range(len(days)):
        day_array = []
        vis = [ name for name in os.listdir(path1) if os.path.isdir(os.path.join(path1, name)) and 'fhd_' in name ]
        vis.sort()
        if visrange != 'all':
            vis = vis[visrange[0]:visrange[1]]
        nocal = 0
        i=-1
        times = []
        cal_dict = {}
        mask_dict = {}
        init_data=False
        for v in vis[0:5]:
            i = i + 1
            obs = v.split('_')[1]
            print('Reading data for obs ' + str(obs))
            data = True
            path = path1 + '/' + v
            if os.path.exists(path + '/calibration/') is False:
                    print('!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!')
                    print('There is no calibration data for this observation. The gains for this observation will be entered as 0.')
                    print(path + '/calibration/')
                    nocal = nocal + 1
                    data = False
                    times.append(None)
            if data == True:
                    obsfile = path + '/metadata/' + obs + '_obs.sav'
                    calfile = path + '/calibration/' + obs + '_cal.sav'
                    settingsfile = path + '/metadata/' + obs + '_settings.txt'
                    cal.read_fhd_cal(calfile,obsfile,settings_file=settingsfile)
            for j in range(nants):
                ant = ant_nums[j]
                if write_calfits == True and data == True:
                    writepath = path + str(nums[i]) + '/calibration/' + obs + '_cal.fits'
                    cal.write_calfits(writepath, clobber=True)
                if i == 0 and data==True:
                    nobs = len(vis)
                    nfreqs = cal.Nfreqs
                    freq_array = cal.freq_array[0,:]
                    freq_array = np.transpose(freq_array)
                    freq_array = np.divide(freq_array,1000000)
                    minfreqind = np.argmax(freq_array>freqrange[0])
                    maxfreqind = np.argmin(freq_array<freqrange[1])
                    cal_dict[ant] = np.empty((nfreqs,nobs))
                    mask_dict[ant] = np.zeros((nfreqs,nobs))
                    if j == 0:
                        ## Get the LST start and end times for this obsid ##
                        time_array = cal.time_array
                        obstime_start = Time(time_array[0],format='jd',location=loc)
                        startTime = obstime_start.sidereal_time('mean').hour
                        JD = int(obstime_start.jd)
                if v == vis[-1] and j == 0:
                    time_array = cal.time_array
                    obstime_end = Time(time_array[-1],format='jd',location=loc)
                    endTime = obstime_end.sidereal_time('mean').hour
                if data is True:
                    if metric == 'amp':
                        gain = np.abs(cal.get_gains(ant, pol))
                    elif metric == 'phase':
                        gain = np.angle(cal.get_gains(ant, pol))
                    else:
                        print('Invalid metric parameter')
                        break
                if data is False:
                    if init_data is False:
                        continue
                    else:
                        gain = np.zeros(np.shape(gain))
                        cal_dict[ant][:,i] = gain[:,0]
                        mask_dict[ant][:,i] = np.ones(np.shape(gain))[:,0]
                        continue
                init_data=True
                cal_dict[ant][:,i] = gain[:,0]
                tArr = cal.time_array
                t = Time(tArr[0],format='jd',location=loc)
                times.append(t)
        mx_dict = {}
        print('Doing per-antenna calculations')
        for j in range(nants):
            ant = ant_nums[j]
            print(ant)
            cal_dict[ant] = np.transpose(cal_dict[ant])
            mask_dict[ant] = np.transpose(mask_dict[ant])
            cal_dict[ant] = cal_dict[ant][:,minfreqind:maxfreqind]
            mask_dict[ant] = mask_dict[ant][:,minfreqind:maxfreqind]
            ## Getting max and min values ##
            calsort = np.sort(cal_dict[ant])
            maxind = int(len(calsort[0])/15)
            maxuse = calsort[:,-maxind]
            maxuse = np.max(maxuse)
            if maxuse >=5:
                    maxind = int(len(calsort[0])/12)
                    maxuse = calsort[:,-maxind]
                    maxuse = np.max(maxuse)
                    print('!!!!!!!!!!!!!!!!!!! Had to adjust max used value !!!!!!!!!!!!!!!!!!!')
            maxval = np.max(calsort[:,-1])
            mx_dict[ant] = np.ma.masked_array(cal_dict[ant], mask=mask_dict[ant])
            freq_array = freq_array[minfreqind:maxfreqind]
            lst_array = get_LSTs(times)
            zenith_RAs = get_zenithRA(times,loc)
            res = {
                    "ant_num": ant,
                    "obsid": obs,
                    "cal_array": cal_dict[ant],
                    "time_array": times,
                    "freq_array": freq_array,
                    "masked_data": mx_dict[ant],
                    "zenith_RA_array": zenith_RAs,
                    "lst_array": lst_array,
                    "pol": pol
                }
            day_array.append(res)
        print('###############################################################')
        print('Out of ' + str(len(vis)) + ' observations, ' + str(nocal) + ' did not contain calibration solutions')
        print('###############################################################')
        data_array[days[k]] = day_array
    return data_array

def getClippedFrequency(freq_array,mx,freq_range=[120,180]):
    idx0 = (np.abs(freq_array - freq_range[0])).argmin()
    idx1 = (np.abs(freq_array - freq_range[1])).argmin()
    freq_array = freq_array[idx0:idx1]
    mx = mx[:,idx0:idx1]
    return mx, freq_array

def getClippedTime(time_array,mx,lst_range=[2,8]):
    idx0 = (np.abs(time_array - np.asarray(lst_range[0]))).argmin()
    idx1 = (np.abs(time_array - np.asarray(lst_range[1]))).argmin()
    time_array = time_array[idx0:idx1]
    mx = mx[idx0:idx1,:]
    return mx, time_array, idx0, idx1