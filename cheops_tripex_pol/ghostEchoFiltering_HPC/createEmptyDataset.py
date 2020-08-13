#! /usr/bin/env python
# coding: utf-8

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import glob as glob
from sys import argv

def emptyDS(timeStart,timeFreq,timeTol,rangeRef,band):
    timeEnd = timeStart + pd.offsets.Hour(1) + pd.offsets.Minute(59) + pd.offsets.Second(59)
    timeEndFirst = timeStart + pd.offsets.Minute(2)
    timeRef = pd.date_range(timeStart, timeEndFirst, freq=timeFreq)
    timeRegrid =  pd.date_range(timeStart, timeEnd, freq=timeFreq)
    refDataset = xr.open_dataset('/scratch2/lterzi/tripex_pol/output/2018/12/01/20181201_00_'+band+'_regridded.nc')
    refDoppler = refDataset.doppler.values
    RadarConst = xr.DataArray(np.empty(len(timeRef))*np.NaN,
                             dims=('time'),
                             coords={'time':timeRef})
    npw1 = xr.DataArray(np.empty(len(timeRef))*np.NaN,
                             dims=('time'),
                             coords={'time':timeRef})
    SNRCorFaCo = xr.DataArray(np.empty([len(timeRef),len(rangeRef)])*np.NaN,
                             dims=('time','range'),
                             coords={'time':timeRef,'range':rangeRef})

    HSDco = xr.DataArray(np.empty([len(timeRef),len(rangeRef)])*np.NaN,
                             dims=('time','range'),
                             coords={'time':timeRef,'range':rangeRef})
     
    SPCco = xr.DataArray(np.empty([len(timeRef),len(rangeRef),len(refDoppler)])*np.NaN,
                             dims=('time','range','doppler'),
                             coords={'time':timeRef,'range':rangeRef,'doppler':refDoppler})
    if band=='ka':
        SNRCorFaCx = xr.DataArray(np.empty([len(timeRef),len(rangeRef)])*np.NaN,
                             dims=('time','range'),
                             coords={'time':timeRef,'range':rangeRef})
        HSDcx = xr.DataArray(np.empty([len(timeRef),len(rangeRef)])*np.NaN,
                             dims=('time','range'),
                              coords={'time':timeRef,'range':rangeRef})
        SPCcx = xr.DataArray(np.empty([len(timeRef),len(rangeRef),len(refDoppler)])*np.NaN,
                             dims=('time','range','doppler'),
                             coords={'time':timeRef,'range':rangeRef,'doppler':refDoppler})


        DataFirst = xr.Dataset({'RadarConst':RadarConst,
                             'npw1':npw1,
                             'SNRCorFaCo':SNRCorFaCo,
                             'SNRCorFaCx':SNRCorFaCx,
                             'HSDco':HSDco,
                             'HSDcx':HSDcx,
                             'SPCco':SPCco,
                             'SPCcx':SPCcx})
        varsName = ['SNRCorFaCo','SNRCorFaCx','HSDco','HSDcx','SPCco','SPCcx']
        DataFirst.time.values[-1]=timeEnd
        DataRegridded = DataFirst[['RadarConst','npw1']].resample(time=timeFreq).nearest(tolerance=timeTol)
        for ind,var in enumerate(varsName):
            DatanewTime = DataFirst[var].resample(time=timeFreq).nearest(tolerance=timeTol)
            if ind == 0:
                Data = xr.merge([DataRegridded,DatanewTime])
            else:
                Data = xr.merge([Data,DatanewTime])

    else:
        DataFirst = xr.Dataset({'RadarConst':RadarConst,
                             'npw1':npw1,
                             'SNRCorFaCo':SNRCorFaCo,
                             'HSDco':HSDco,
                             'SPCco':SPCco})
        varsName = ['SNRCorFaCo','HSDco','SPCco']
    
        DataFirst.time.values[-1]=timeEnd
        DataRegridded = DataFirst[['RadarConst','npw1']].resample(time=timeFreq).nearest(tolerance=timeTol)    
        for ind,var in enumerate(varsName):
            if (var == 'SPCco'):
                 dataEmpty = np.empty((len(timeRegrid),len(rangeRef),len(refDoppler)))
                 DatanewTime = xr.DataArray(dataEmpty,
                             dims=('time','range','doppler'),
                             coords={'time':timeRegrid,'range':rangeRef, 'doppler':refDoppler})
                 DatanewTime.name = var
                 for tind,ti in enumerate(timeRegrid):
                     DataStep = DatanewTime[tind,:,:]*np.NaN
                     DatanewTime[tind,:,:] = DataStep
            else:
                DatanewTime = DataFirst[var].resample(time=timeFreq).nearest(tolerance=timeTol)
            if ind == 0:
                Data = xr.merge([DataRegridded,DatanewTime])
            else:
                Data = xr.merge([Data,DatanewTime])


    return Data
