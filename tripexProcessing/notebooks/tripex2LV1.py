import glob
import numpy as np
import pandas as pd
import xarray as xr
from sys import argv
import glob
from netCDF4 import Dataset

import writeData
import tripexLib as trLib
import readRadarInfo as rdInfo

#-- code to bring data to LVL1 as described in Dias neto et al. 2019 --> this code is supposed to regrid the X,Ka,W-band radars to the same time and height grid

#-- completely new version, because Joses skript needs a lot of input variables:
# variables provided by controll file:
scriptname, inputPathX,inputPathKa,inputPathW, dayIndex, prefixL1, outputPath, timeFreq, timeTolerance, beginRangeRef, endRangeRef, rangeFreq, rangeTolerance, rangeGateOffset, zeOffset = argv

# define range grid:
beginRangeRef = int(beginRangeRef)
endRangeRef = int(endRangeRef)
rangeFreq = int(rangeFreq)
rangeTolerance = int(rangeTolerance)
rangeRef = np.arange(beginRangeRef, endRangeRef, rangeFreq)
usedIndexRange = np.ones((len(rangeRef)))*np.nan
rangeGateOffSet = float(rangeGateOffset)
zeOffSet = float(zeOffset)

# select correct date (on cheops we want to process all days at the same time, it is easier to just provide an index)
file_date = open("datelist.txt")
dateList = file_data.readlines()
dayIndex = int(dayIndex)
date2proc = dateList[dayIndex]

# define the time grid:
#timeFreq = 'S'
start = (' ').join([date2proc,'00:00:00'])
end = (' ').join([date2proc,'23:59:59'])
timeRef = pd.date_range(start,end, freq=timeFreq)

#-- define the outputfile
#prefixL1 = 'tripex_pol_triple_freq_LVL1'
outputFile = ('_').join([prefixL1, dateName, date2proc+'.nc'])
outPutFilePath = ('/').join([outputPath, outputFile])

#- read in and regrid the data to a common time height grid:
filePathX = '/'.join([inputPathX,date2proc+'*.nc'])
filePathKa = '/'.join([inputPathKa,date2proc+'*.nc'])
filePathW = '/'.join([inputPathW,date2proc+'*.nc'])

filesX = glob.glob(filePathX)
filesKa = glob.glob(filePathKa)
filesW = glob.glob(filePathW)

for fX,fKa,fW in zip(filesX,filesKa,filesW):
    dataX = xr.open_dataset(fX)
    dataXnewTime = dataX.resample(time=timeRef).nearest(tolerance=timeTolerance)
    dataXRegridded = dataXnewTime.reindex(range=newRange,method='nearest',tolerance=rangeTolerance)
    dataKa = xr.open_dataset(fKa)
    dataKanewTime = dataKa.resample(time=timeRef).nearest(tolerance=timeTolerance)
    dataKaRegridded = dataKanewTime.reindex(range=newRange,method='nearest',tolerance=rangeTolerance)
    dataW = xr.open_dataset(fW)
    dataWnewTime = dataW.resample(time=timeRef).nearest(tolerance=timeTolerance)
    dataWRegridded = dataWnewTime.reindex(range=newRange,method='nearest',tolerance=rangeTolerance)
    
    #-- add offset correction to Ka Band
    #-- check if Ze is linear or dB
    #-- write netCDF file with new data




