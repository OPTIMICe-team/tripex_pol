#####################################################
# This script is meant to regrid the polarimetric Moments that were measured during the tripex-pol Campaign. It also recalculates KDP, because the automatically calculated one is wrong.
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

dateStart = pd.to_datetime('20190122'); dateEnd = pd.to_datetime('20190221')
dateList = pd.date_range(dateStart, dateEnd,freq='D')
dataOutPath = '/work/lvonterz/tripex_pol/output/'

for dayIndex,date2proc in enumerate(dateList[0:-1]):
    print('now reprocessing '+date2proc.strftime('%Y%m%d'))
    date2end = dateList[dayIndex+1]
    dataPath = '/data/obs/campaigns/tripex-pol/wband_scan/l0/Y{year}/M{month}/D{day}'.format(year=date2proc.strftime('%Y'),month=date2proc.strftime('%m'),day=date2proc.strftime('%d'))
    fileId = '*{dateStr}_??????_P09_CEL.LV1.nc'.format(dateStr=date2proc.strftime('%y%m%d'))
    filePath = os.path.join(dataPath, fileId)
    fileList = sorted(glob.glob(filePath))

    #-- merge all files into one dataset:
    oneData = pro.createOneDataset(fileList)
    #-- concatenate all chirps into one range coordinate
    mergedData = pro.mergeChirps(oneData)
    #-- now I need to recalulate PhiDP and KDP:
    PhiDPData = pro.calcPhiDP(mergedData)
    KDPData = pro.calcKDP(PhiDPData)

    #-- define reference grid (so rangeRef and timeRef) according to non-pol data:
    timeFreq = '4S'; timeTol = '2S';
    date2end = date2end - pd.offsets.Second(4)
    timeRef = pd.date_range(date2proc,date2end,freq=timeFreq)
    beginRangeRef=0; endRangeRef=12000; rangeFreq=36; rangeTol=18
    rangeRef = np.arange(beginRangeRef,endRangeRef,rangeFreq)
    dataRegridded = pro.regridPol(KDPData,rangeRef,timeRef,rangeTol,timeTol)
    #-- now add the correct attributes:
    dataRegridded = pro.globalAttr(dataRegridded)
    dataRegridded = pro.varAttr(dataRegridded)    
    #-- now write the data to netcdf:
    dataOutName = '{datestr}_tripex_pol_poldata_L1_mom.nc'.format(datestr=date2proc.strftime('%Y%m%d'))
    dataRegridded.to_netcdf(dataOutPath+dataOutName)
    #plotOutPath = '/work/lvonterz/tripex_pol/plotting/plots/pol/'
    #plotID = 'LVL1_pol_mom'
    #plot.plotPol(dataRegridded,plotOutPath,date2proc.strftime('%Y%m%d'),plotID)
