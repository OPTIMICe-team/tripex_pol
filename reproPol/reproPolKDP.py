#####################################################
# This script is meant to regrid the KDP that was calculated from Alexander. We dont have data for all the hours...
#####################################################

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import glob as glob
from sys import argv
from netCDF4 import Dataset

#own routines:
import processData as pro
#import plotting_routines as plot

def createOneDataset(fileList):
    finalData = xr.Dataset()
    for filePath in fileList:
        data = xr.open_dataset(filePath)
        data.Time.attrs['units']='Seconds since 01.01.2001 00:00:00'
        data = xr.decode_cf(data)
        # get rid of duplicate time values:
        _, index_time = np.unique(data['Time'], return_index=True)
        data = data.isel(Time=index_time)
        finalData = xr.merge([finalData,data])
    finalData = finalData.rename({'KDP':'KDP_Alexander','Range':'range'})
    return finalData

def createEmptyDataArray(timeRef,rangeRef):
    data = xr.DataArray(np.empty([len(timeRef),len(rangeRef)])*np.nan,
                        dims = ('time','height'),
                        coords = {'time':timeRef,'height':rangeRef})
    data = data.rename('KDP_Alexander')
    return data

####################################################################################################################
dateStart = pd.to_datetime('20181101'); dateEnd = pd.to_datetime('20190221')
dateList = pd.date_range(dateStart, dateEnd,freq='D')
dataOutPath = '/work/lvonterz/tripex_pol/output/'
dataLV1InPath = '/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_1/'
for dayIndex,date2proc in enumerate(dateList[0:-1]):
    print('now reprocessing '+date2proc.strftime('%Y%m%d'))
    date2end = dateList[dayIndex+1]
    dataPath = '/work/lvonterz/tripex_pol/input/CEL/{year}/{month}/{day}'.format(year=date2proc.strftime('%Y'),month=date2proc.strftime('%m'),day=date2proc.strftime('%d'))
    fileId = '*{dateStr}_??????_P09_CEL.LV1.nc'.format(dateStr=date2proc.strftime('%y%m%d'))
    filePath = os.path.join(dataPath, fileId)
    fileList = sorted(glob.glob(filePath))
    if fileList:# only if we have files on that day
        #-- merge all files into one dataset:
        oneData = createOneDataset(fileList)
        #-- define reference grid (so rangeRef and timeRef) according to non-pol data:
        timeFreq = '4S'; timeTol = '2S';
        date2end = date2end - pd.offsets.Second(4)
        timeRef = pd.date_range(date2proc,date2end,freq=timeFreq)
        beginRangeRef=0; endRangeRef=12000; rangeFreq=36; rangeTol=18
        rangeRef = np.arange(beginRangeRef,endRangeRef,rangeFreq)
        dataRegridded = pro.regridPol(oneData,rangeRef,timeRef,rangeTol,timeTol)
        #-- now add the correct attributes:
        dataRegridded = pro.varAttr(dataRegridded,KDP=True)
        dataRegridded = dataRegridded.transpose('time','height')
    else:# if we dont have files on that day, add empty dataset to LV1 file
        print('no KDP data for that day, now creating empty dataset')
        timeFreq = '4S'; timeTol = '2S';
        date2end = date2end - pd.offsets.Second(4)
        timeRef = pd.date_range(date2proc,date2end,freq=timeFreq)
        beginRangeRef=0; endRangeRef=12000; rangeFreq=36; rangeTol=18
        rangeRef = np.arange(beginRangeRef,endRangeRef,rangeFreq)
        dataRegridded = createEmptyDataArray(timeRef,rangeRef)
        dataRegridded.attrs['standard_name'] = 'KDP'
        dataRegridded.attrs['long_name'] = 'Specific differential phase shift'
        dataRegridded.attrs['units'] = 'deg km-1'
        dataRegridded = dataRegridded.transpose('time','height')
    #-- now add KDP to the LV1 dataset:
    dataLV1Name = '{datestr}_tripex_pol_poldata_L1_mom.nc'.format(datestr=date2proc.strftime('%Y%m%d'))
    if glob.glob(dataLV1InPath+dataLV1Name):
        dataLV1 = xr.open_dataset(dataLV1InPath+dataLV1Name)
        dataLV1 = xr.merge([dataLV1,dataRegridded])
        dataLV1.to_netcdf(dataOutPath+dataLV1Name)



