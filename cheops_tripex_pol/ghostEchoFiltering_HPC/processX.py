#############################################################################################
#this skript removes masks the X-band spectra according to the fast and slow edges of the Ka-Band
#############################################################################################

import os
import numpy as np
import pandas as pd
import xarray as xr
import glob as gb
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import scipy

def MaskXr(xData, kaData, maxVelKa, minVelKa, kaSpec, kadoppler, channel, xSpec=0, xdoppler=0):
# the funtion name is a bit weird, but basically this just transfers the max, min velocities as well as the x-spec into xarray datasets to make the handling easier later on
    maxVelKaXR= xr.DataArray(maxVelKa, 
                             coords={'range':kaData.range},
                             dims=('range'))
    minVelKaXR= xr.DataArray(minVelKa, 
                             coords={'range':kaData.range}, 
                             dims=('range'))
    
    time = kaData.time.values
    time = time.reshape(1)
    SpecKa = xr.DataArray(kaSpec, 
                          dims=('time','range','doppler'),
                          coords={'time':time,'range':kaData.range, 'doppler':kadoppler})
    #plt.plot(SpecKa.doppler,SpecKa[0,25,:])
    #plt.show()
    if channel == 'o':
        SpecX = xr.DataArray(xSpec, 
                          dims=('time','range','doppler'), 
                          coords={'time':time,'range':xData.range, 'doppler':xdoppler})
        return maxVelKaXR, minVelKaXR, SpecX, SpecKa
    else:
        return maxVelKaXR, minVelKaXR, SpecKa
    
def spectraCleaning(xSpec, maxVelKa, minVelKa):
# now this is where the magic happens: the X-Band spectra are masked according to the Ka-band fast and slow edges. I have to do this in a loop over the ranged unfortunately, because it would require too much memory if removed directly
    datasets = []
    for rInd,rangeSel in enumerate(maxVelKa.range):
        maxVel = maxVelKa[rInd].values
        minVel = minVelKa[rInd].values
        DataFiltered = xSpec.loc[dict(range=rangeSel)].where(xSpec.doppler < maxVel) 
        DataFiltered = DataFiltered.where(xSpec.doppler >= minVel)
        datasets.append(DataFiltered)
    maskedSpec = xr.concat(datasets, dim='range')
    maskedSpec = maskedSpec.transpose('time','range','doppler')
    return maskedSpec

