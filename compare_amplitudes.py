#!/usr/bin/python

from astropy.io import fits
from pyuvdata import UVData
import numpy as np
import healpy as hp
import sys
import math
import matplotlib.pyplot as plt
import scipy.io

#### SET PARAMETERS ####

dirty_origin = 'FHD' #Allows for reading in dirty visibilities from a separate uvfits file (assuming it would be coming from pyuvsim)
data_directory = 'SingleOffZenithSourceTest'
factor2normalization = False #Scale up the model visibilities by a factor of 2 - useful when comparing against pyuvsim

path = '/Users/home/Documents/_Files/Dara/School/Graduate/RadCos/FHD_Pyuvsim_comparison/AllGaussian/' + data_directory + '/'
name = 'ref_1.1_offzenith_gauss'
prefix = path + name + '_'
simpath = '/Users/home/Documents/_Files/Dara/School/Graduate/RadCos/FHD_Pyuvsim_comparison/' + name + '.uvfits'
pol = 'XX' #Which polarization you want to look at

autoRenorm = False #Try automatically calculating and applying a renormalization factor between the two visibility sets
##########################


files = [prefix + f for f in ['flags.sav','vis_XX.sav','params.sav', 'vis_YY.sav','settings.txt','vis_model_XX.sav','vis_model_YY.sav']]

dirty = UVData()
model = UVData()

#Model is the fhd simulation, dirty is the pyuvsim simulation
model.read_fhd(files, use_model=True)

if dirty_origin == 'pyuvsim':
    dirty.read(simpath)
    path = path + 'PYVIS_'
elif dirty_origin == 'FHD':
    dirty.read_fhd(files, use_model=False)
else:
    print('ERROR: Unknown pyuvsim origin')
path = path + pol + '_'

mdat = model.get_data(pol)
mdat = np.abs(mdat)
if factor2normalization == True:
    path = path + 'factor2renorm_'
    mdat = np.multiply(mdat, 2)
ddat = dirty.get_data(pol)
ddat = np.abs(ddat)

#Plot a histogram of the ratio between visibilities
ratio = np.divide(mdat, ddat)
print(ratio)
ratplot = plt.figure()
plt.hist(ratio,30,[0,2])
plt.xlabel('FHD/Pyuvsim Visibility Amplitude (' + pol + ' polarization)')
plt.ylabel('Count')
plt.title('Histogram of Ratio of Visibility Amplitudes')
plt.show()
ratplot.savefig(path + 'ratioPlot.pdf')

test = np.histogram(ratio,30)

ind = np.argmax(test[0])
#a = test[1][ind]
#b = test[1][ind+1]

#Get an estimate of the ratio between the normalizations of the dirty and model visibilities
mult = np.average([test[1][ind],test[1][ind+1]])
print(mult)

separateplot = plt.figure()
plt.hist(ddat,50,[0,7],alpha=0.5,label='pyuvsim')
plt.hist(mdat,50,[0,7],alpha=0.5,label='fhd')
plt.legend(loc='upper right')
plt.xlabel('Visibility Amplitude (' + pol + ')')
plt.ylabel('Count')
plt.title('Histogram of Visibility Amplitudes')
plt.show()
separateplot.savefig(path + 'separatePlot.pdf')

sep_fullrange = plt.figure()
#maxes = [np.maximum(ddat),np.maximum(mdat)]
#print(maxes)
max = np.maximum(np.amax(ddat),np.amax(mdat))
#print(max)
plt.hist(ddat,50,[0,max],alpha=0.5,label='pyuvsim')
plt.hist(mdat,50,[0,max],alpha=0.5,label='fhd')
plt.legend(loc='upper right')
plt.xlabel('Visibility Amplitude (' + pol + ')')
plt.ylabel('Count')
plt.title('Histogram of Visibility Amplitudes')
plt.show()
sep_fullrange.savefig(path + 'sep_fullrange.pdf')

if autoRenorm == True:
    ddatadj = np.multiply(ddat, mult)
    sep_postAdjustment = plt.figure()
    plt.hist(ddatadj,50,[0,max],alpha=0.5,label='pyuvsim')
    plt.hist(mdat,50,[0,max],alpha=0.5,label='fhd')
    plt.legend(loc='upper right')
    plt.xlabel('Visibility Amplitude (' + pol + ')')
    plt.ylabel('Count')
    plt.title('Histogram of Visibility Amplitudes After Rescaling')
    plt.show()
    sep_postAdjustment.savefig(path + 'sep_postAdjustment.pdf')

## Try adjusting based on max value instead of center of distribution
if autoRenorm == True:
    topmult = np.divide(np.amax(mdat),np.amax(ddat))
    ddatadj = np.multiply(ddat, topmult)
    adjustByMax = plt.figure()
    plt.hist(ddatadj,50,[0,max],alpha=0.5,label='pyuvsim')
    plt.hist(mdat,50,[0,max],alpha=0.5,label='fhd')
    plt.legend(loc='upper right')
    plt.xlabel('Visibility Amplitude (' + pol + ')')
    plt.ylabel('Count')
    plt.title('Histogram of Visibility Amplitudes After Rescaling By Max Value')
    plt.show()
    adjustByMax.savefig(path + 'adjustByMax.pdf')

#Adjust the visibilities so that their normalizations are equivalent now
#dirty.data_array = np.multiply(dirty.data_array,mult)
#dirty.write_uvfits('/Users/home/Documents/_Files/Dara/School/Graduate/RadCos/FHD_Pyuvsim_comparison/Results/ref_1.1_uniform_pyuvsim_adjustedAmplitude.uvfits', spoof_nonessential=True)
