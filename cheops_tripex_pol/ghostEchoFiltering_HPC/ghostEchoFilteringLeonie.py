#! /usr/bin/env python
# coding: utf-8
#######################################################################################
# this is the main script for the ghostEcho Filtering. This manages all parts necessary to generate the LV0 and LV1 files for the tripex-pol campaign
# 1) regridding X,Ka and W Band data to same grid
# 2) removing ghostEchos and saving all spectra to one ncfile
# 3) calculating the LV1 files (so the moments)
#######################################################################################

import os
import time 
import numpy as np
import pandas as pd
import xarray as xr
import glob as glob
from sys import argv

#import multiprocessing as mp
#import matplotlib.pyplot as plt
from netCDF4 import Dataset
import readData as rd
import process as pro
import kaLib
import wLib
import regridGhostecho as regrid
import processX as proX
import moments as mo
import createEmptyDataset as createDataset
########################################################################################
#-- the dayindex that is supposed to be processed is input from the bash skript
scriptname, dayIndex = argv
dayIndex = int(dayIndex)

################################## 1) now we are reading in the data and regridding everything #################################
#-- define range reference:
beginRangeRef=0; endRangeRef=12000; rangeFreq=36; rangeTolerance=18 
rangeOffsetKa = 2.2; rangeOffsetX = 0.32+36 # !!!we also have to add 36 until 14th of december to the X-Band!!!
#-- because our W-band has the chirps: we need to define 3 reference heights (I would have loved to merge the chirps before doing the regridding, but that is extremely memory consuming)
#I need a different rangeRef for the days where we only have 3 chirps (so some days in november, I dont remember which ones that were,. definitely 01.11.)
beginRangeC2 = 1476; beginRangeC3 = 3996; beginRangeC4 = 8172  #beginRangeC2 = 1224; beginRangeC3 = 5508 
#-- create reference range vectors
rangeRef = np.arange(beginRangeRef,endRangeRef,rangeFreq)
rangeRefC1 = np.arange(beginRangeRef,beginRangeC2,rangeFreq)
rangeRefC2 = np.arange(beginRangeC2,beginRangeC3,rangeFreq)
rangeRefC3 = np.arange(beginRangeC3,beginRangeC4,rangeFreq)#,beginRangeC4,rangeFreq)
rangeRefC4 = np.arange(beginRangeC4,endRangeRef,rangeFreq)
#-- create reference time vectors
timeFreq = '4S'; timeTolerance = '2S'; 
dateStart = pd.to_datetime('20181101'); dateEnd = pd.to_datetime('20190221')
dateList = pd.date_range(dateStart, dateEnd,freq='D')
date2proc = dateList[dayIndex]; date2end = dateList[dayIndex+1]
timeRef = pd.date_range(date2proc,date2end,freq=timeFreq)
#-- define output path
Output = '/scratch2/lterzi/tripex_pol/output/'
OutputPath = os.path.join(Output,
                            date2proc.strftime('%Y'),
                            date2proc.strftime('%m'),
                            date2proc.strftime('%d'))

#-- regrid all data to the same time-range grid with afore defined rangeRef and timeRef:
regrid.regridOneDay(date2proc,OutputPath,timeFreq,rangeRef,rangeRefC1,rangeRefC2,rangeRefC3,rangeRefC4,timeTolerance,rangeTolerance,rangeOffsetKa,rangeOffsetX)


############################################# 2) now the ghostEcho Filtering part ###################################
#-- filter X and W band data according to spectra of Ka band:
print(' ')
print('Now starting with ghostEchoFiltering')
print('')

kapath = sorted(glob.glob(OutputPath+'/*_ka_regridded.nc'))
xpath = sorted(glob.glob(OutputPath+'/*_x_regridded.nc'))
wpath = sorted(glob.glob(OutputPath+'/*_w_regridded.nc'))

ind = 0
fileIDdate = date2proc.strftime('%Y%m%d')
file_last = OutputPath+'/'+fileIDdate+'_23_LVL0_data_regridded.nc'
files_already_there = sorted(glob.glob(OutputPath+'/*_LVL0_data_regridded.nc'))

#because kadata and xdata were saved in 2 hourly files, and wdata in 1 hourly files, I am now reading in every wdata file and every second loop a new ka and x data file
if not file_last in files_already_there:
    if any(files_already_there):# now if the processing was killed somewhere in the middle, I am finding out the last file that was processed. Only works when the filename is as it is on Cheops for now...
        file_last_processed = files_already_there[-1]
        HourProcessed = float(file_last_processed[-30:-28])
        HourProcessed = HourProcessed + 1.0 
        Hours2proc = np.arange(HourProcessed,24.0)
    else:
        HourProcessed = 0.0
        Hours2proc = np.arange(0.0,24.0)

    # now process the hours that were not processed before:
    for hour in Hours2proc:
        hour2proc = '{:02.0f}'.format(hour)
        print(hour2proc)
        #-- now read in the regridded w-band file: 
        wDataFile = OutputPath+'/'+fileIDdate+'_'+hour2proc+'_w_regridded.nc'
        if wDataFile in wpath:
            wData = xr.open_dataset(wDataFile)
        else: # if we dont have data for this hour, I will create an empty Dataset for that hour: (I should put that as a funtion in emptyDataset.py, but it is working now and I dont wanna change anything ;) )
            timeStart = date2proc + pd.offsets.Hour(int(hour2proc))
            timeEnd = timeStart + pd.offsets.Minute(59) + pd.offsets.Second(59)
            timeRef = pd.date_range(timeStart, timeEnd, freq=timeFreq)
            DoppLen = np.asarray([512, 512, 512, 512])
            MaxVel = np.asarray([10.2643, 10.2643, 10.2643, 10.2643])
            
            C1HNoisePow = xr.DataArray(np.empty([len(timeRef),len(rangeRefC1)]),
                                       dims=('Time','C1Range'),
                                       coords={'Time':timeRef,'C1Range':rangeRefC1})
            C2HNoisePow = xr.DataArray(np.empty([len(timeRef),len(rangeRefC2)]),
                                       dims=('Time','C2Range'),
                                       coords={'Time':timeRef,'C2Range':rangeRefC2})
            C3HNoisePow = xr.DataArray(np.empty([len(timeRef),len(rangeRefC3)]),
                                       dims=('Time','C3Range'),
                                       coords={'Time':timeRef,'C3Range':rangeRefC3})
            C4HNoisePow = xr.DataArray(np.empty([len(timeRef),len(rangeRefC4)]),
                                       dims=('Time','C4Range'),
                                       coords={'Time':timeRef,'C4Range':rangeRefC4})
            
            C1VNoisePow = xr.DataArray(np.empty([len(timeRef),len(rangeRefC1)]),
                                       dims=('Time','C1Range'),
                                       coords={'Time':timeRef,'C1Range':rangeRefC1})
            C2VNoisePow = xr.DataArray(np.empty([len(timeRef),len(rangeRefC2)]),
                                       dims=('Time','C2Range'),
                                       coords={'Time':timeRef,'C2Range':rangeRefC2})
            C3VNoisePow = xr.DataArray(np.empty([len(timeRef),len(rangeRefC3)]),
                                       dims=('Time','C3Range'),
                                       coords={'Time':timeRef,'C3Range':rangeRefC3})
            C4VNoisePow = xr.DataArray(np.empty([len(timeRef),len(rangeRefC4)]),
                                       dims=('Time','C4Range'),
                                       coords={'Time':timeRef,'C4Range':rangeRefC4})
            
            
            C1HSpec = xr.DataArray(np.empty([len(timeRef),len(rangeRefC1),DoppLen[0]]),
                                   dims=('Time','C1Range','C1Vel'),
                                   coords={'Time':timeRef,'C1Range':rangeRefC1})
            C2HSpec = xr.DataArray(np.empty([len(timeRef),len(rangeRefC2),DoppLen[0]]),
                                   dims=('Time','C2Range','C2Vel'),
                                   coords={'Time':timeRef,'C2Range':rangeRefC2})
            C3HSpec = xr.DataArray(np.empty([len(timeRef),len(rangeRefC3),DoppLen[0]]),
                                   dims=('Time','C3Range','C3Vel'),
                                   coords={'Time':timeRef,'C3Range':rangeRefC3})
            C4HSpec = xr.DataArray(np.empty([len(timeRef),len(rangeRefC4),DoppLen[0]]),
                                   dims=('Time','C4Range','C4Vel'),
                                   coords={'Time':timeRef,'C4Range':rangeRefC4})
            C1VSpec = xr.DataArray(np.empty([len(timeRef),len(rangeRefC1),DoppLen[0]]),
                                   dims=('Time','C1Range','C1Vel'),
                                   coords={'Time':timeRef,'C1Range':rangeRefC1})
            C2VSpec = xr.DataArray(np.empty([len(timeRef),len(rangeRefC2),DoppLen[0]]),
                                   dims=('Time','C2Range','C2Vel'),
                                   coords={'Time':timeRef,'C2Range':rangeRefC2})
            C3VSpec = xr.DataArray(np.empty([len(timeRef),len(rangeRefC3),DoppLen[0]]),
                                   dims=('Time','C3Range','C3Vel'),
                                   coords={'Time':timeRef,'C3Range':rangeRefC3})
            C4VSpec = xr.DataArray(np.empty([len(timeRef),len(rangeRefC4),DoppLen[0]]),
                                   dims=('Time','C4Range','C4Vel'),
                                   coords={'Time':timeRef,'C4Range':rangeRefC4})
            DoppLen = xr.DataArray(DoppLen)
            MaxVel = xr.DataArray(MaxVel)
            wData = xr.Dataset({'C1HNoisePow':C1HNoisePow,
                                'C2HNoisePow':C2HNoisePow,
                                'C3HNoisePow':C3HNoisePow,
                                'C4HNoisePow':C4HNoisePow,
                                'C1VNoisePow':C1VNoisePow,
                                'C2VNoisePow':C2VNoisePow,
                                'C3VNoisePow':C3VNoisePow,
                                'C4VNoisePow':C4VNoisePow,
                                'C1HSpec':C1HSpec,
                                'C2HSpec':C2HSpec,
                                'C3HSpec':C3HSpec,
                                'C4HSpec':C4HSpec,
                                'C1VSpec':C1VSpec,
                                'C2VSpec':C2VSpec,
                                'C3VSpec':C3VSpec,
                                'C4VSpec':C4VSpec,
                                'DoppLen':DoppLen,
                                'MaxVel':MaxVel,
                                })
    
        # I dont know why I did it in this complicated way...
        DataStart = pd.to_datetime(wData.Time.values[0])
        DataHour = DataStart.strftime('%H')
        DataHourFloat = float(DataHour)
        if (DataHourFloat == 0) or ((DataHourFloat % 2) == 0): # now every second hour read in ka and x band
            kaDataFile = OutputPath+'/'+fileIDdate+'_'+hour2proc+'_ka_regridded.nc'
            xDataFile = OutputPath+'/'+fileIDdate+'_'+hour2proc+'_x_regridded.nc'
            if kaDataFile in kapath:
                kaData = xr.open_dataset(kaDataFile)
            else:# if no file exists, create empty dataset:
                timeStart = date2proc + pd.offsets.Hour(int(hour2proc))
                kaData = createDataset.emptyDS(timeStart,timeFreq,timeTolerance,rangeRef,'ka')
            if xDataFile in xpath:
                xData = xr.open_dataset(xDataFile)
            else:
                timeStart = date2proc + pd.offsets.Hour(int(hour2proc))
                xData = createDataset.emptyDS(timeStart,timeFreq,timeTolerance,rangeRef,'x') 
        print('kaData',kaData)
        print('xData',xData)
        print('wData',wData)
      
#-- now finally the Ghostecho filtering, the Data is then later saved to a netCDF file (noise and Spec)
        regrid.ghostEchoFiltering(kaData,xData,wData,OutputPath,rangeRef)



###################################################################
#3) now calculate the moments from the regridded and filtered data:

fileID = '_LVL0_data_regridded_test'
fileIDdate = date2proc.strftime('%Y%m%d')
dataPath = OutputPath+'/*'+fileID+'.nc'
dataOutPath = os.path.join(OutputPath, fileIDdate+'_moments.nc')

#-- calculating the moments 
fileList = sorted(glob.glob(dataPath))
momentsW = xr.Dataset()
momentsX = xr.Dataset()
momentsKa = xr.Dataset()

for fileName in fileList:
# calculate the moments and merge them into one file per day:
    data = xr.open_dataset(fileName)
    tempMomeW = mo.getMoments(data,'W')
    momentsW = xr.merge([momentsW,tempMomeW])
    tempMomeKa = mo.getMoments(data,'Ka')
    momentsKa = xr.merge([momentsKa,tempMomeKa])
    tempMomeX = mo.getMoments(data,'X')
    momentsX = xr.merge([momentsX,tempMomeX])

#-- global attributes     
momentsX = mo.globalAttributes(momentsX)
momentsX = mo.variableAtrributes(momentsX,'X')  
momentsKa = mo.globalAttributes(momentsKa)
momentsKa = mo.variableAtrributes(momentsKa,'Ka')
momentsW = mo.globalAttributes(momentsW)
momentsW = mo.variableAtrributes(momentsW,'W')

moments = xr.merge([momentsX,momentsKa,momentsW])
#encDic = {'XZeH':{'zlib':True},
#          'XZeV':{'zlib':True},
#          'WZeH':{'zlib':True},
#          'WZeV':{'zlib':True},
#          'KaZeH':{'zlib':True},
#          'KaZeV':{'zlib':True},
#          'Xrv':{'zlib':True},
#          'Karv':{'zlib':True},
#          'Wrv':{'zlib':True},
#          'sk':{'zlib':True},
#          'sw':{'zlib':True},
#          'time':{'zlib':True},
#          'range':{'zlib':True}}

moments.to_netcdf(dataOutPath)
