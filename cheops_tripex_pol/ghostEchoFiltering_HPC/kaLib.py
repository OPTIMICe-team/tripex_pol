################################################################################
# this script was written by Jose, I made slight adjustments, but mainly performance wise. This script contains all the functions used to calculate Ka-X-spectra and noise, aswell as the edges of the Ka-spectra
###############################################################################

import os
import numpy as np
import pandas as pd
import xarray as xr
import glob as gb
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def getCalSpect(dataSet, ch):
# function from Jose, I added notes to that which might not be correct...    
#ch = o or x (co channel or cross channel)
    
    rangeLen = dataSet.range.shape[0]
    try:
        timeLen = dataSet.time.shape[0]
    except:
        timeLen = 1
    
    rangeData = dataSet.range.values.reshape(1, rangeLen, 1)
    #-- radar Constant related to 5km Height and 200ns pulses (description from nc-file)
    radConst = dataSet.RadarConst.values.reshape(timeLen,1,1)
    #-- SNRCorFa: Signal-to-noise-factor corrected based on Receiver calibration
    SNRCorFa = dataSet['SNRCorFaC'+ch].values.reshape(timeLen,rangeLen,1)
    #-- Doppler Spectrum
    SPC = dataSet['SPCc'+ch].values
    #-- Noise power
    npw = dataSet.npw1.values.reshape(timeLen,1,1)
    #-- provided by Stefan    
    calSPC = radConst*SNRCorFa*(rangeData**2/5000.**2)*SPC/npw
    
    return calSPC


def getCalNoise(dataSet, ch):
    
    rangeLen = dataSet.range.shape[0]
    try:
        timeLen = dataSet.time.shape[0]
    except:
        timeLen = 1
    
    rangeData = dataSet.range.values.reshape(1, rangeLen, 1)
    radConst = dataSet.RadarConst.values.reshape(timeLen,1,1)
    SNRCorFa = dataSet['SNRCorFaC'+ch].values.reshape(timeLen,rangeLen,1)
    #-- ???? What does HSD mean????
    HSD = dataSet['HSDc'+ch].values.reshape(timeLen, rangeLen, 1)
    npw = dataSet.npw1.values.reshape(timeLen,1,1)
    
    calNoise = radConst*SNRCorFa*(rangeData**2/5000.**2)*HSD/npw
    
    return calNoise


def reorderSpec(dopplerVal, specData, dpMin, dpMax):
# Joses skript!    
#-- reorder Spectrum around maximum value (-> switch left and right side of spectrum)
    dpIndex = np.arange(len(dopplerVal))
    
    maxDpIndex = dpIndex[dopplerVal >= 0].max()
    zeroIndex = dpIndex[dopplerVal >= 0].min()
    
    newDoppler = np.ones_like(dopplerVal)
    newSpec = np.ones_like(specData) * np.nan
    
    newDoppler[len(newDoppler) - maxDpIndex-1:] = dopplerVal[:maxDpIndex+1]
    newDoppler[:len(newDoppler) - maxDpIndex-1] = dopplerVal[maxDpIndex+1:]
    
    newSpec[:,:,len(newDoppler) - maxDpIndex-1:] = specData[:,:,:maxDpIndex+1]
    newSpec[:,:,:len(newDoppler) - maxDpIndex-1] = specData[:,:,maxDpIndex+1:]

    newDoppler = newDoppler*(-1)
    try:
        doppLimIndMin = dpIndex[newDoppler>=dpMax].max()
    except:
        doppLimIndMin = 0
    try:
        doppLimIndMax = dpIndex[newDoppler<=dpMin].min()
    except:
        doppLimIndMax = len(dopplerVal)
    
    newDoppler = newDoppler[doppLimIndMin:doppLimIndMax+1]
    newSpec = newSpec[:,:,doppLimIndMin:doppLimIndMax+1]
   
    return newDoppler, newSpec


def getMaxAndMinSpecVel(spectra, velArr):
# Joses script!
#-- calculates max and min velocity    
    velMatrix = (spectra/spectra)*velArr
    maxVel = np.nanmax(velMatrix, axis=1)
    minVel = np.nanmin(velMatrix, axis=1)

    return maxVel, minVel


def loadKaData(dataPathKa, wTime):
#- Joses script!
#- opens dataset from Ka band at that specific hour
    if pd.to_datetime(wTime).hour % 2 == 0:

        fileId = pd.to_datetime(wTime).strftime('%Y%m%d_%H')
    
    else:
    
        fileId = pd.to_datetime(wTime) - pd.to_timedelta('1H')
        fileId = fileId.strftime('%Y%m%d_%H')
    
    filePathKa = os.path.join(dataPathKa, fileId+'*.znc')
    fileListW = gb.glob(filePathKa)

    kaData = xr.open_dataset(fileListW[0])
    kaData.time.values = kaData.time.values + kaData.microsec.values*10**(-6)
    kaData.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
    kaData = xr.decode_cf(kaData)
    
    return kaData


def removSuspData(calNoise, newSpecKa, noiseCutoff):
#- joses script, use edge velocities at 3db above noise?

    noiseLevelArr = 10*np.log10(calNoise[0,:,0])+noiseCutoff

    for index, noiseLevel in enumerate(noiseLevelArr):

        specRow = newSpecKa[0][index]
        specRow[specRow < noiseLevel] = np.nan
        newSpecKa[0][index] = specRow
        
    return newSpecKa


def getMaskProf(kaData, newSpecKa, newDopplerKa, minHeight, Nyquist):
#- Joses script!
#- calculate max and min velocities    
    maxVelKa, minVelKa = getMaxAndMinSpecVel(newSpecKa[0], newDopplerKa)

    zeros = np.zeros_like(kaData.range.values)
    maxVelKa = np.nansum(np.array([zeros,maxVelKa]), axis=0)
    minVelKa = np.nansum(np.array([zeros,minVelKa]), axis=0)

    maxVelKa[kaData.range.values<minHeight] = np.nan
    minVelKa[kaData.range.values<minHeight] = np.nan
    maxVelKa[maxVelKa>= Nyquist] = Nyquist + 10
    minVelKa[minVelKa<=-Nyquist] = Nyquist - 10

    return maxVelKa, minVelKa


