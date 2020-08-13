import glob
import numpy as np
import pandas as pd

import xarray as xr
from sys import argv
import glob
import os
from netCDF4 import Dataset
import time

#-- import own functions:
import readData as rd

outputPath = '/work/lvonterz/tripexPro/outputtest/'
inputPathKa = '/data/obs/site/jue/joyrad35/'
inputPathW = '/data/obs/campaigns/tripex-pol/wband_gra/l0/'
inputPathX = '/data/obs/site/jue/joyrad10/'
date2proc = '20181208'
date2end = '20181209'
prefix = 'LVL1_Test_Ka'
#- define reference time grid
timeFreq = '4S'
timeTolerance = '2S'
date = pd.to_datetime(date2proc)
#- define reference range grid
beginRangeRef=100
endRangeRef=12000
rangeFreq=35
rangeTolerance=15
varVecKa = ['SNRCorFaCo','SPCco']
varVecW = ['MaxVel','DoppLen','Elv']
varVecX = ['Doppler','SPCco','npw1','SNRCorFaCo','HSDCo']
 #'C1VNoisePow''C2VNoisePow','C3VNoisePow','C4VNoisePow','C1VSpec','C2VSpec', 'C3VSpec','C4VSpec'-- they are not passed to the function, because they depend on the chirp number... therefore they are directly included in the function
rangeRef = np.arange(beginRangeRef, endRangeRef, rangeFreq)
timeRef = pd.date_range(date2proc,date2end,freq=timeFreq)

ZeOffset = 5.0
#-- Ka input file:

filePathKa = os.path.join(inputPathKa,
                         date.strftime('%Y'),
                         date.strftime('%m'),
                         date.strftime('%d'))
filesKa = glob.glob(filePathKa+'/*.znc')
for indfile,fKa in enumerate(filesKa[0:1]):
    kaDataRegridded = rd.regridKaXData(fKa,varVecKa,timeFreq,rangeRef, timeTolerance,rangeTolerance)
    print(kaDataRegridded)
    quit()
'''
filePathW = os.path.join(inputPathW,
                        date.strftime('%Y'),
                         date.strftime('%m'),
                         date.strftime('%d'))
filesW = glob.glob(filePathW+'/*LV0.nc')
for indfile,fW in enumerate(filesW[0:1]):
    wDataRegridded = rd.regridWData(fW,varVecW,timeRef,rangeRef, timeTolerance,rangeTolerance)
    print(wDataRegridded.C1Range)
    print(wDataRegridded.C2Range)

filePathX = os.path.join(inputPathX,
                        date.strftime('%Y'),
                         date.strftime('%m'),
                         date.strftime('%d'))
filesX = glob.glob(filePathX+'/*.znc')
for indfile,fX in enumerate(filesX[0:1]):
    XDataRegridded = rd.regridKaXData(fX,varVecX,timeFreq,rangeRef, timeTolerance,rangeTolerance)
'''    
















