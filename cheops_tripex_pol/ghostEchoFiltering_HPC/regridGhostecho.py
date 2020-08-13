####################################################################################################
# now this script contains both the "managment" for regridding the LV0 data and the ghostEcho Filtering. I dont know why I have them both in here, but it doesnt make a difference I guess.. ;) 
####################################################################################################
import os
import time 
import numpy as np
import pandas as pd
import xarray as xr
import glob as gb
from sys import argv

import multiprocessing as mp
#import matplotlib.pyplot as plt
from netCDF4 import Dataset
import readData as rd
import process as pro
import kaLib
import wLib
import processX as proX

def regridOneDay(date,OutputPath,timeRef,rangeRef,
                rangeRefC1,rangeRefC2,rangeRefC3,rangeRefC4,
                timeTol,rangeTol,rangeOffsetKa,rangeOffsetX):
# this function reads in and regrids the datasets and then saves it to a netcdf file. the function calls the extern function readData.py

    files_already_there = gb.glob(OutputPath+'/*.nc')
    iDDate = date.strftime('%Y%m%d')
    file_last_w = OutputPath+'/'+iDDate+'_23_w_regridded.nc'
    file_last_ka = OutputPath+'/'+iDDate+'_22_ka_regridded.nc'
    file_last_x = OutputPath+'/'+iDDate+'_22_x_regridded.nc'
     
    print('processing {date}'.format(date=date)) 
    inputPathW = '/scratch2/lterzi/tripex_pol/wband/'
    inputPathKa = '/scratch2/lterzi/tripex_pol/kaband/'
    inputPathX = '/scratch2/lterzi/tripex_pol/xband/'
    
    dataPathW = os.path.join(inputPathW, 
                             date.strftime('%Y'), 
                             date.strftime('%m'), 
                             date.strftime('%d'), 
                             '*LV0.nc')

    fileListW = sorted(gb.glob(dataPathW))

    dataPathKa = os.path.join(inputPathKa,
                             date.strftime('%Y'),
                             date.strftime('%m'),
                             date.strftime('%d'))
    fileListKa = sorted(gb.glob(dataPathKa+'/*.znc'))
    dataPathX = os.path.join(inputPathX,
                             date.strftime('%Y'),
                             date.strftime('%m'),
                             date.strftime('%d'))
    fileListX = sorted(gb.glob(dataPathX+'/*.znc'))
    
    varVecW = ['MaxVel','DoppLen','Elv']
    varVecX = ['SPCco','SNRCorFaCo','HSDco']
    varVecKa = ['SNRCorFaCo','SNRCorFaCx','HSDco','HSDcx','SPCco','SPCcx']

    #-- first of all regrid all data (X,Ka,W) to the same grid
     
    for indFileKa,filePathKa in enumerate(fileListKa):
        print(filePathKa)
        kaDataRegridded, DataHour = rd.regridKaData(date,filePathKa,varVecKa,timeRef,rangeRef,timeTol,rangeTol,rangeOffsetKa)
        kaDataRegridded.to_netcdf(OutputPath+'/'+iDDate+'_'+DataHour+'_ka_regridded.nc')

    
    for indFileX,filePathX in enumerate(fileListX):
        print(filePathX)
        xDataRegridded, DataHour = rd.regridXData(date,filePathX,varVecX,timeRef,rangeRef,timeTol,rangeTol,rangeOffsetX)
        xDataRegridded.to_netcdf(OutputPath+'/'+iDDate+'_'+DataHour+'_x_regridded.nc')
    
    if file_last_w not in files_already_there:
        for indFileW,filePathW in enumerate(fileListW):
            wDataRegridded, DataHour = rd.regridWData(date, filePathW,varVecW,timeRef,rangeRefC1,rangeRefC2,rangeRefC3,rangeRefC4,timeTol,rangeTol)
            wDataRegridded.to_netcdf(OutputPath+'/'+iDDate+'_'+DataHour+'_w_regridded.nc')
        
       
def ghostEchoFiltering(kaData2proc,xData2proc,wData2proc,outputPath,rangeRef):
# now this function makes all the ghostecho filtering happening. First I am merging the chirps from w-band into one range coordinate
    for wIndex,wTime in enumerate(wData2proc.Time.values):
        kaData = kaData2proc.sel(time=wTime)
        xData = xData2proc.sel(time=wTime)
        wDataOld = wData2proc.sel(Time=wTime)
        maxVel = wDataOld.MaxVel.values
        doppLen = wDataOld.DoppLen.values
        dv_vec = np.empty(4)
        for chirp in range(4):
            ChRange = 'C{chirp}Range'.format(chirp=chirp+1)
            ChSpec = 'C{chirp}HSpec'.format(chirp=chirp+1)
            ChSpecV = 'C{chirp}VSpec'.format(chirp=chirp+1)
            ChVel = 'C{chirp}Vel'.format(chirp=chirp+1)
            chNoiseH = 'C{chirp}HNoisePow'.format(chirp=chirp+1)
            chNoiseV = 'C{chirp}VNoisePow'.format(chirp=chirp+1)

            velData = np.linspace(-maxVel[chirp], maxVel[chirp], doppLen[chirp], dtype=np.float32)
            dv_diff = np.diff(velData)
            dv_vec[chirp] = dv_diff[0]
            wDataOld = wDataOld.assign({ChVel:velData})
            velRef = np.linspace(-maxVel[0], maxVel[0], doppLen[0], dtype=np.float32) 
            wDataOld = wDataOld.reindex({ChVel:velRef}, method = 'nearest', tolerance = 0.05)
            wDataOld[ChSpecV] = wDataOld[ChSpecV]/dv_vec[chirp]
            wDataOld[ChSpec] = wDataOld[ChSpec]/dv_vec[chirp]
            if chirp == 0:
                DataTime1 = wDataOld[[ChSpecV,ChSpec]]#, chNoiseH, chNoiseV]]
                DataTime1 = DataTime1.rename({ChVel:'Vel',ChRange:'range',ChSpecV:'SpecV',ChSpec:'SpecH'})#,chNoiseH:'NoiseH',chNoiseV:'NoiseV'})
                wData = DataTime1
            
            else:
                DataTime1 = wDataOld[[ChSpecV,ChSpec]]#,chNoiseH, chNoiseV]]
                DataTime1 = DataTime1.rename({ChVel:'Vel',ChRange:'range',ChSpecV:'SpecV',ChSpec:'SpecH'})#,chNoiseH:'NoiseH',chNoiseV:'NoiseV'})
                wData = xr.concat([wData,DataTime1],dim='range')
        print(wData)
       

        # now most of the functions used here are directly used from Jose
        for channel in ['o','x']:
            #-- calculate the right ka spectrum and noise:
            calSpect = kaLib.getCalSpect(kaData, channel)
            calNoise = kaLib.getCalNoise(kaData, channel)
            specCalNoiCorrKa = calSpect - (np.ones_like(calSpect)*calNoise)
            #-- the ka spectrum is stored as 0,..,10,-10,..0, so we need to reorder the spectrum 
            newDopplerKa, newSpecKa = kaLib.reorderSpec(kaData.doppler.values,
                                                        specCalNoiCorrKa,
                                                        -10.56824, 
                                                        10.56824)# NyquistVelocity
            

            newSpecKa = 10*np.log10(newSpecKa)
            newSpecKa = kaLib.removSuspData(calNoise, 
                                            newSpecKa, 3)
            # get the maximum and minimum velocity to then mask w and x band
            maxVelKa, minVelKa = kaLib.getMaskProf(kaData, newSpecKa, 
                                                   newDopplerKa, rangeRef[0],
                                                   10.56824)
             
            if channel == 'x':
                # X-band doesnt have x-channel, so we dont need to calculate it here!
                maxVelKaXr, minVelKaXr, xrSpecKa = proX.MaskXr(xData, kaData, maxVelKa, minVelKa, newSpecKa,newDopplerKa, channel) 
                #-- now mask the W-spectrum according to fast and slow edge from ka-band
                maskedSpecW = pro.spectraCleaning(wData, maxVelKaXr, minVelKaXr,channel)
                #- calculate Ze, MDV
                dv=dv_vec[0]
                ZeKa = xrSpecKa.sum(dim='doppler')
                ZeW = maskedSpecW['SpecH'].sum(dim='Vel')*dv
                mdvKa = ((xrSpecKa*(xrSpecKa['doppler'])).sum(dim='doppler'))/ZeKa
                mdvW = ((maskedSpecW['SpecH']*(maskedSpecW['Vel'])).sum(dim='Vel')*dv)/ZeW
                #-- write output to nc file
                pro.writeOutPut(wTime, wIndex, outputPath, maskedSpecW, calNoise, xrSpecKa, channel, ZeKa, mdvKa, ZeW, mdvW,dv) 
                print('Done with spectra cleaning, timestep', maskedSpecW.Time, 'now writing output')
            else:
                # now also x-band, so first calculate it and remove noise and reorder it:
                calSpectX = kaLib.getCalSpect(xData, channel)
                calNoiseX = kaLib.getCalNoise(xData, channel)
                specCalNoiCorrX = calSpectX - (np.ones_like(calSpectX)*calNoiseX)
                newDopplerX, newSpecX = kaLib.reorderSpec(xData.doppler.values,
                                                        specCalNoiCorrX,
                                                        -20.0, 
                                                        20.0)# NyquistVelocity
                newSpecX = 10*np.log10(newSpecX)
                newSpecX = kaLib.removSuspData(calNoiseX, 
                                                newSpecX, 3)
                #- mask the X- and W-band:
                maxVelKaXr, minVelKaXr, xrSpecX, xrSpecKa = proX.MaskXr(xData, kaData,
                                                     maxVelKa, minVelKa, newSpecKa,newDopplerKa, channel, newSpecX, newDopplerX) 
                maskedSpecX = proX.spectraCleaning(xrSpecX,maxVelKaXr,minVelKaXr)
                maskedSpecW = pro.spectraCleaning(wData, maxVelKaXr, minVelKaXr,channel)
                #- calculate Ze, MDV
                dv = dv_vec[0]
                ZeKa = (10**(xrSpecKa/10)).sum(dim='doppler')
                ZeW = maskedSpecW['SpecV'].sum(dim='Vel')*dv
                ZeX = (10**(maskedSpecX/10)).sum(dim='doppler')
                
                mdvKa = ((xrSpecKa*(xrSpecKa['doppler'])).sum(dim='doppler'))/ZeKa
                mdvW = ((maskedSpecW['SpecV']*(maskedSpecW['Vel'])).sum(dim='Vel')*dv)/ZeW
                mdvX = ((maskedSpecX*(maskedSpecX['doppler'])).sum(dim='doppler'))/ZeX

         
                print('Done with spectra cleaning, timestep', maskedSpecW.Time, 'now writing output')
                pro.writeOutPut(wTime, wIndex, outputPath, maskedSpecW, calNoise,  xrSpecKa, channel, ZeKa, mdvKa, ZeW, mdvW,dv, ZeX, mdvX, calNoiseX, maskedSpecX) 

    




















