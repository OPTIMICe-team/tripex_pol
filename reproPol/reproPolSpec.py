#####################################################
# This script is meant to regrid the polarimetric spectra that were measured during the tripex-pol Campaign. 
#####################################################

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import glob as glob
from sys import argv
from netCDF4 import Dataset
from matplotlib import pyplot as plt
#own routines:
import processData as pro
#import plotting_routines as plot
#############################################################################
def getNewNipySpectral():
    
    from matplotlib import cm
    from matplotlib.colors import ListedColormap

    numEnt = 15
    
    viridis = cm.get_cmap('nipy_spectral', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    
    colorSpace = np.linspace(198, 144, numEnt)/256
    colorTest=np.zeros((numEnt,4))
    colorTest[:,3] = 1
    colorTest[:,0]=colorSpace

    newcolors[- numEnt:, :] = colorTest
    newcmp = ListedColormap(newcolors)
    
    return newcmp
############################################################################

debugging = True

dateStart = pd.to_datetime('20190122'); dateEnd = pd.to_datetime('20190131')
dateList = pd.date_range(dateStart, dateEnd,freq='D')
dataOutPath = '/work/lvonterz/tripex_pol/output/'

for dayIndex,date2proc in enumerate(dateList):
    date2end = dateList[dayIndex+1]-pd.offsets.Second(4)
    date2proc = date2proc + pd.offsets.Hour(15)
    timeFreq = '1H'; timeHour = pd.date_range(date2proc,date2end,freq=timeFreq)
    dataPath = '/data/obs/campaigns/tripex-pol/wband_scan/l0/Y{year}/M{month}/D{day}'.format(year=date2proc.strftime('%Y'),
                                                                                             month=date2proc.strftime('%m'),
                                                                                             day=date2proc.strftime('%d'))
    for hourInd,hour in enumerate(timeHour):
        #the files are stored according to the hour-minute-second they begun. So when one set started the hour before, but still goes on in this hour, we need to also include the data from the hour before, so that this is not lost during the regridding: (maybe there is a better way, but for now it workes and is not too slow)
        h1 = hour - pd.offsets.Hour(1)
        dataPath1 = '/data/obs/campaigns/tripex-pol/wband_scan/l0/Y{year}/M{month}/D{day}'.format(year=h1.strftime('%Y'),
                                                                                             month=h1.strftime('%m'),
                                                                                             day=h1.strftime('%d'))
 
        fileId1 = '*{dateStr}_{HH}????_P09_CEL.LV0.nc'.format(dateStr=h1.strftime('%y%m%d'),HH=h1.strftime('%H'))
        filePath = sorted(glob.glob(os.path.join(dataPath1, fileId1)))
        if filePath:
            fileList = [filePath[-1]]
        else:
            fileList = []
        fileId2 = '*{dateStr}_{HH}????_P09_CEL.LV0.nc'.format(dateStr=date2proc.strftime('%y%m%d'),HH=hour.strftime('%H'))
        filePath2 = os.path.join(dataPath, fileId2)
        #depending on the memory available on your computer, you may need to change this to fileList=glob.glob(filePath2) if we have a continuous dataset (so no scans were made in the hour and the data is available during the entire hour)
        fileList.extend(sorted(glob.glob(filePath2)))
        if fileList: #- only process if files are there
            print('now processing: ',hour)
            if debugging==True:
                print(fileList)
            #-- merge all files of this hour into one dataset:
            data = pro.createOneDataset(fileList,spec=True)
            
            
            #-- now merge the chirps into one range, and while we are at it, I am going to recalculate the VSpec (noise doesnt need to be removed, this is already done automatically in he RPG software) 
            # (the VSpec that is in this dataset for now is not really the VSpec, but a composite of H and V, and with the ChReVHSpec variable we can calculate the real VSpec. 
            if debugging==True:
                print('now merging chirps')
            dataMerged = pro.mergeChirps(data,spec=True)
            dataMerged['sZDR'] = 10*np.log10(dataMerged['HSpec']) - 10*np.log10(dataMerged['VSpec'])
            dataMerged['sSNR_H'] = 10*np.log10(dataMerged['sSNR_H'])
            dataMerged['sSNR_V'] = 10*np.log10(dataMerged['sSNR_V'])
       

            #-- make smallest velocities of spectrum 0 (because we are looking at an angle, the doppler velocity is influenced by the horizontal velocity. As we dont really now how large that is, and since Alexander is also just setting Vel_min to 0, I am going to do the same
            if debugging==True:
                print('calculating max and min vel')
            maxVelH,minVelH = pro.calcOffset(dataMerged,'HSpec')
            dataMerged['maxVelH'] = maxVelH
            dataMerged['minVelH'] = minVelH
            maxVelV,minVelV = pro.calcOffset(dataMerged,'VSpec')
            dataMerged['maxVelV'] = maxVelV
            dataMerged['minVelV'] = minVelV
            maxVelZDR,minVelZDR = pro.calcOffset(dataMerged,'sZDR')
            dataMerged['maxVelZDR'] = maxVelZDR
            dataMerged['minVelZDR'] = minVelZDR
            #-- now subtract maxvel values from vel array, that should put right edge to 0:
            if debugging==True:
                print('calculating offset and removing it')
            dataOffset = pro.removeOffset(dataMerged)
            

            #-- if you want to use the data at this stage, you should mask the data with SNR_H > 10, but for the regridding I don't need it. There is going to be a comment that the data should not be used if SNR_H < 10.0! Same goes for the spectral data (so sSNR)
            #dataMasked = dataOffset.where(10*np.log10(dataOffset['SNR_H']) > 10.0)
            #dataMasked = dataOffset.where(dataOffset['sSNR_H'] > 10.0)


            #-- now, lets regrid this dataset to the same grid as the non-pol dataset
            timeFreq = '4S'; timeTol = '2S'; time2end = hour + pd.offsets.Hour(1) - pd.offsets.Second(4)
            timeRef = pd.date_range(hour, time2end,freq=timeFreq)
            beginRangeRef=0; endRangeRef=12000; rangeFreq=36; rangeTol=18
            rangeRef = np.arange(beginRangeRef,endRangeRef,rangeFreq)
            if debugging==True:
                print('regrid Dataset')
            dataRegridded = pro.regridPol(dataOffset,rangeRef,timeRef,rangeTol,timeTol,spec=True)
            #dataRegridded['ZDR_ma'] = dataRegridded['ZDR'].where(10*np.log10(dataRegridded['SNR_H'])>10)


            #-- now lets write the attributes for our dataset 
            data2save = pro.varAttr(dataRegridded,spec=True)
            data2save = pro.globalAttr(data2save)
            #-- and save the dataset
            if debugging==True:
               print('save the file')
            dataOutName = '{datestr}_tripex_pol_poldata_L0_spec_regridded.nc'.format(datestr=hour.strftime('%Y%m%d_%H'))
            dataOutputPath = dataOutPath + '{year}/{month}/{day}/'.format(year=date2proc.strftime('%Y'),
                                                                          month=date2proc.strftime('%m'),
                                                                          day=date2proc.strftime('%d')) 
            # we need this so that the ncfile is not too large (if we leave this out, than the ncfile is 7.5GB per hour!!!!)
            encDic = {'HSpec':{'zlib':True},
                      'VSpec':{'zlib':True},
                      'HNoisePow':{'zlib':True},
                      'VNoisePow':{'zlib':True},
                      'SNR_H':{'zlib':True},
                      'SNR_V':{'zlib':True},
                      'sSNR_H':{'zlib':True},
                      'sSNR_V':{'zlib':True},
                      'ZH':{'zlib':True},
                      'ZV':{'zlib':True},
                      'ZDR':{'zlib':True},
                      'sZDR':{'zlib':True},
                      'maxVelH':{'zlib':True},
                      'maxVelV':{'zlib':True},
                      'maxVelZDR':{'zlib':True},
                      'Vel2ZeroH':{'zlib':True},
                      #'Vel2ZeroV':{'zlib':True},
                      #'Vel2ZeroZDR':{'zlib':True},
                      'time':{'zlib':True},
                      'height':{'zlib':True}}
            data2save.to_netcdf(dataOutputPath+dataOutName,encoding=encDic)
        else:
            print('no files available for ',hour)
    #quit()


