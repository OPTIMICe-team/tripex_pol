#################################################################################################
# this skript contains all the functions needed to reprocess the polarimetric dataset. For now this calculated KDP from RHO_HV like Jose did it, but as Alexander said, this is not the best way. 
# this script also conatins functions to move the fastest edge of the spectra to 0
# author: Leonie von Terzi
#################################################################################################

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import glob as glob
from sys import argv
from netCDF4 import Dataset
import matplotlib.pyplot as plt


def createOneDataset(fileList,spec=False,KDP=False):
# this functions simply merges all ncfiles from one hour into one dataset, so that we have a simpler format to later regrid the data
# input: name of files that need to be merged, and if you set spec=True, it merges the spectral data (I am doing it like this right now, because I don't need all the variables thT re in the ncfile, and if I select them specifically I save a lot of memory and disk space
    finalData = xr.Dataset()
    for filePath in fileList:
        data = xr.open_dataset(filePath)
        # the time stored in the ncfiles is not the real time, but we have the milliseconds stored extra, so I need to manually decode it:
        data.Time.values = data.Time.values + data.Timems.values/1000.
        data.Time.attrs['units']='Seconds since 01.01.2001 00:00:00'
        data = xr.decode_cf(data)
        # get rid of duplicate time values:
        _, index_time = np.unique(data['Time'], return_index=True)
        data = data.isel(Time=index_time)
        chirpNum = data.ChirpNum.values
        if spec==True: # for the polarimetric spectra
            for chirp in range(chirpNum):
                chirpName = 'C{chirp}'.format(chirp=chirp+1)
                #deltaRange = data[chirpName+'Range']/data[chirpName+'Range']*data.RangeRes[chirp]
                temp = xr.Dataset({chirpName+'Range':data[chirpName+'Range'],
                                   chirpName+'HSpec':data[chirpName+'HSpec'],
                                   chirpName+'VSpec':data[chirpName+'VSpec'],
                                   chirpName+'ReVHSpec':data[chirpName+'ReVHSpec'],
                                   chirpName+'HNoisePow':(data[chirpName+'HNoisePow']),
                                   chirpName+'VNoisePow':data[chirpName+'VNoisePow'],
                                   })
                finalData = xr.merge([finalData, temp])
            temp = xr.Dataset({'Azm':data.Azm,
                               'Elv':data.Elv,
                               'ChirpNum':data.ChirpNum,
                               'RangeRes':data.RangeRes,
                               'MaxVel':data.MaxVel,
                               'DoppLen':data.DoppLen,
                               })
        else:
            for chirp in range(chirpNum):# for the LVL1 polarimetric moments:
                chirpName = 'C{chirp}'.format(chirp=chirp+1)
                #deltaRange = data[chirpName+'Range']/data[chirpName+'Range']*data.RangeRes[chirp]
                temp = xr.Dataset({chirpName+'Range':data[chirpName+'Range'],
                                   chirpName+'ZDR':data[chirpName+'ZDR'],
                                   chirpName+'PhiDP':data[chirpName+'PhiDP'],
                                   chirpName+'RHV':data[chirpName+'RHV'],
                                   chirpName+'ZE':10*np.log10(data[chirpName+'ZE']),
                                   chirpName+'Skew':data[chirpName+'Skew'],
                                   chirpName+'SpecWidth':data[chirpName+'SpecWidth'],
                                   chirpName+'MeanVel':data[chirpName+'MeanVel'],
                                   chirpName+'Kurt':data[chirpName+'Kurt'],
                                   })
                finalData = xr.merge([finalData, temp])

            temp = xr.Dataset({'Azm':data.Azm,
                               'Elv':data.Elv,
                               'ChirpNum':data.ChirpNum,
                               'RangeRes':data.RangeRes})

        finalData = xr.merge([finalData, temp])
    return finalData

def mergeChirps(data,spec=False):
# this function merges the chirps of the dataset. I am doing this, because the non-polarimetric X and Ka-Band don't have a chirps and therefore we only have one fixed range-resolution
# input: data as xarray-dataset, if you set spec=True, the spectral data will be merged.
    chirpNum = data.ChirpNum.values
    if spec==True:
        maxVel = data.MaxVel.values; doppLen = data.DoppLen.values; dv_vec = np.empty(chirpNum)
        for chirp in range(chirpNum):
            ChRange = 'C{chirp}Range'.format(chirp=chirp+1)
            ChHSpec = 'C{chirp}HSpec'.format(chirp=chirp+1)
            ChVSpec = 'C{chirp}VSpec'.format(chirp=chirp+1)
            ChReVHSpec = 'C{chirp}ReVHSpec'.format(chirp=chirp+1)
            ChHNoisePow = 'C{chirp}HNoisePow'.format(chirp=chirp+1)
            ChVNoisePow = 'C{chirp}VNoisePow'.format(chirp=chirp+1)
            ChVel = 'C{chirp}Vel'.format(chirp=chirp+1)
            
            #- calculate noise density:
            NoiseDensV = data[ChVNoisePow]/doppLen[chirp]
            NoiseDensH = data[ChHNoisePow]/doppLen[chirp]
            #- now we need to decompose the VSpec, because this is not actually the vert. Spectrum but saved as a composite of H and V
            data[ChVSpec] = 2*data[ChVSpec] - data[ChHSpec] - 2*data[ChReVHSpec]
            #- there was a software mistake, so the noise of ReHV was not stored correctly. Therefore Alexander suggested to use s SNR threshold of 10dB, otherwise the data will be masked. For this we need to calculate SNR: SNR = signal power/noise power. In order to calculate the correct values, we need to mask all values below -90 dBZ.
            specThreshold = 10**(-90/10)
            data[ChHSpec] = data[ChHSpec].where(data[ChHSpec]>specThreshold,np.NaN)
            data[ChVSpec] = data[ChVSpec].where(data[ChVSpec]>specThreshold,np.NaN)
            NoisePowH = data[ChHSpec].count(dim=ChVel)*NoiseDensH
            NoisePowV = data[ChVSpec].count(dim=ChVel)*NoiseDensV
            SignalPowH = data[ChHSpec].sum(dim=ChVel)
            SignalPowV = data[ChVSpec].sum(dim=ChVel)
            data['SNR_H'] = SignalPowH / NoisePowH
            data['SNR_V'] = SignalPowV / NoisePowV
            #-- now also make it spectral (so spectral ZDR against SNR in respective bin (we see also some high ZDR values at the left (fast) edge, which should not be there)) and then the data should be masked with sSNR > 10dB aswell
            data['sSNR_H'] = data[ChHSpec] / NoiseDensH
            data['sSNR_V'] = data[ChVSpec] / NoiseDensH 
            #-- it is the easiest to calculate Z,ZDR,ZDP at this stage, because once the chirps are merged, it is difficult to sort out the dv needed to integrate over the spectrum...
            data['ZH'] = data[ChHSpec].sum(dim=ChVel)
            data['ZV'] = data[ChVSpec].sum(dim=ChVel)
            data['ZDR'] = 10*np.log10(data['ZH'])-10*np.log10(data['ZV'])
            data['ZDP'] = data['ZH']-data['ZV']
            #- because the different chirps have a different Doppler resolution, we also need to regrid along that axis:
            velData = np.linspace(-maxVel[chirp], maxVel[chirp], doppLen[chirp], dtype=np.float32) # in the dataformat, the dopplervelocity was not assigned yet, but rather it is stored as maxVel and doppLen which you then need to manually assign to the doppler Vel coordinate
            dv_diff = np.diff(velData) # since we regrid along the doppler axis, we need to divide the regridded doppler spectra by dv
            dv_vec[chirp] = dv_diff[0]
            data = data.assign({ChVel:velData})
            velRef = np.linspace(-maxVel[0], maxVel[0], doppLen[0], dtype=np.float32)# just use the Doppler velocity from the smallest chirp
            data = data.reindex({ChVel:velRef}, method = 'nearest', tolerance = 0.05) # regrid
            data[ChVSpec] = data[ChVSpec]/dv_vec[chirp]
            data[ChHSpec] = data[ChHSpec]/dv_vec[chirp]
            data[ChReVHSpec] = data[ChReVHSpec]/dv_vec[chirp]
            #- now we can rename the variables to without the chirps and then we merge the datasets along the range coordinate
            dataCh = data[[ChVSpec,ChHSpec,ChHNoisePow,ChVNoisePow,'ZH','ZV','ZDR','ZDP','SNR_H','SNR_V','sSNR_H','sSNR_V']]
            dataCh = dataCh.rename({ChRange:'range',ChVSpec:'VSpec',ChHSpec:'HSpec',ChHNoisePow:'HNoisePow',ChVNoisePow:'VNoisePow',ChVel:'Vel'})
            if chirp==0:
                finalData = dataCh
            else:
                finalData = xr.concat([finalData,dataCh],dim='range')
        delRange = np.concatenate((selfDiv(data.C1Range)*data.RangeRes[0],
                           selfDiv(data.C2Range)*data.RangeRes[1],
                           selfDiv(data.C3Range)*data.RangeRes[2],))
        delRange = xr.DataArray(delRange,
                            dims=('range'),
                            coords={'range':finalData.range})
        dv =  xr.DataArray(dv_vec,dims=('Chirp'))
        temp = xr.Dataset({'Azm':data.Azm,
                       'Elv':data.Elv,
                       'RangeRes':data.RangeRes,
                       'delRange':delRange,
                       'dv':dv})
    else:
        for chirp in range(chirpNum):
            ChRange = 'C{chirp}Range'.format(chirp=chirp+1)
            ChZDR = 'C{chirp}ZDR'.format(chirp=chirp+1)
            ChPhiDP = 'C{chirp}PhiDP'.format(chirp=chirp+1)
            ChRHV = 'C{chirp}RHV'.format(chirp=chirp+1)
            ChZe = 'C{chirp}ZE'.format(chirp=chirp+1)
            ChSkew = 'C{chirp}Skew'.format(chirp=chirp+1)
            ChWidth = 'C{chirp}SpecWidth'.format(chirp=chirp+1)
            ChKurt = 'C{chirp}Kurt'.format(chirp=chirp+1)
            ChVel = 'C{chirp}MeanVel'.format(chirp=chirp+1)
            # now merge the dataset along the range coordinate
            dataCh = data[[ChZDR,ChPhiDP,ChRHV,ChZe,ChSkew,ChWidth,ChKurt,ChVel]]
            dataCh = dataCh.rename({ChRange:'range',ChZDR:'ZDR',ChPhiDP:'PhiDP',ChRHV:'RHV',ChZe:'DBZ',ChSkew:'SK',ChWidth:'WIDTH',ChKurt:'KURT',ChVel:'VEL'})
            if chirp==0:
                finalData = dataCh    
            else:
                finalData = xr.concat([finalData,dataCh],dim='range')
    
        delRange = np.concatenate((selfDiv(data.C1Range)*data.RangeRes[0],
                           selfDiv(data.C2Range)*data.RangeRes[1],
                           selfDiv(data.C3Range)*data.RangeRes[2],))
        delRange = xr.DataArray(delRange,
                            dims=('range'),
                            coords={'range':finalData.range})
        temp = xr.Dataset({'Azm':data.Azm,
                       'Elv':data.Elv,
                       'RangeRes':data.RangeRes,
                       'delRange':delRange})
    
    finalData = xr.merge([finalData, temp])
    return finalData 
    
def calcPhiDP(data):
# copied from Jose directly
    PHIDP = data['PhiDP'].values
    for i in range(len(PHIDP)):
        PHIDP[i][~np.isnan(PHIDP[i])] = np.unwrap(PHIDP[i][~np.isnan(PHIDP[i])])
    PHIDP = PHIDP*-1
    PhiDPxr = xr.DataArray(PHIDP, 
                       dims=('Time', 'range'), 
                       coords={'Time':data.Time,
                               'range':data.range})
    PhiDPxr = np.rad2deg(PhiDPxr)
    PhiDPxr = PhiDPxr.rolling(range=5, min_periods=1,center=True).mean()
    PhiDPxr.attrs = {'standard_name':'PhiDP',
                     'long_name': 'Differential phase shift',
                     'units':'deg'}
    PhiDPxr = PhiDPxr.rename('PhiDP')
    data = data.drop('PhiDP')
    data = xr.merge([data,PhiDPxr])
    return data

def selfDiv(var):
    return var/var

def calcKDP(data):
# taken from Jose
    delRange = data.delRange.values
    timeWindow = 60 # so jose has here 7*7, but his timeref  is 7s resolution, mine is 4, so to get same averaginf, I have to multiply by 60
    PHIDP = data['PhiDP'].rolling(Time=timeWindow,min_periods=1,center=True).mean()
    PHIDP = PHIDP.rolling(range=5,min_periods=3, center=True).mean()
    KDP = PHIDP.diff(dim='range')/(2.*abs(delRange[0:-1])*1e-3)
    KDP = KDP.rename('KDP')
    KDP.attrs = {'long_name': 'Specific differential phase shift',
                 'units':'deg/km'}
    data = xr.merge([data,KDP])
    return data

def regridPol(data,rangeRef,timeRef,rangeTol,timeTol):
# this regrids the polarimetric dataset to the same time/height grid as the non-polarimetric dataset:
# input: data as xarray-dataset, range and time reference vectors as well as the tolerances to use.
    #-- since we are looking slanted, we need to project the range coordinate to a "height" coordinate:
    try:
        height = data.range.values*np.sin(np.deg2rad(data.Elv.values[0]))
    except:
        height = data.range.values*np.sin(np.deg2rad(30.0))
    data.range.values = height
    data = data.rename({'range':'height','Time':'time'})
    _, index_time = np.unique(data['height'], return_index=True)
    data = data.isel(height=index_time)

    dataNewRange = data.reindex({'height':rangeRef}, method='nearest',tolerance=rangeTol)  
    dataNewTime = dataNewRange.reindex({'time':timeRef}, method='nearest', tolerance=timeTol)
    dataNewTime.height.attrs = {'long_name':'projected height',
                                'units':'m'}
    return dataNewTime
    
def calcMom(data):
# this calculates the moments from the spectra, but we already do this in the merging of the chirps, so this is just to test the routines..
    ZeH = data['HSpec'].sum(dim='vel')*data.dv
    ZeV = data['VSpec'].sum(dim='vel')*data.dv
    ZDR = 10*np.log10(ZeH)-10*np.log10(ZeV)
    ZDP = ZeH-ZeV
    moments = xr.Dataset({'ZeH':ZeH,'ZeV':ZeV,'ZDR':ZDR,'ZDP':ZDP})
    return moments

def calcOffset(data,var):
# this function calculates the right and left edge of the polarimetric dataset (so meaning that the vector we get here can be used to move the fast edge of the psectra to 0). This is needed, because since we are looking slanted, the Doppler velocities are influenced by the horizontal wind component as well ( so not just the vertical as it would be if we look zenith). Thespectra therefore looks more snake like and we don't have the actual information of the fall velocities of the particles. For now moving it to 0 is the easiest and cleanest way to deal with it. 
    vel = data.Vel.values
    spectra = data[var].values
    velMatrix = spectra/spectra*vel
    maxVel = np.nanmax(velMatrix,axis=2)
    minVel = np.nanmin(velMatrix,axis=2)
    minVelXR = xr.DataArray(minVel,
                            dims=('Time','range'),
                            coords={'Time':data.Time,'range':data.range})
    maxVelXR = xr.DataArray(maxVel,
                            dims=('Time','range'),
                            coords={'Time':data.Time,'range':data.range})
    
    return maxVelXR,minVelXR
  
def removeOffset(data):
# here the offset we calculated in calcOffset is used to retrieve a matrix with which we can move the fast edge to 0 m/s. The dataset needs to hae the variable maxVel from th calcOffset function.
    newVelH = data.Vel - data['maxVelH'].fillna(0)
    newVelH = newVelH.transpose('Time','range','Vel')
    data['Vel2ZeroH'] = newVelH.fillna(0)
#    newVelV = data.Vel - data['maxVelV'].fillna(0)
#    newVelV = newVelV.transpose('Time','range','Vel')
#    data['Vel2ZeroV'] = newVelV
#    newVelZDR = data.Vel - data['maxVelZDR'].fillna(0)
#    newVelZDR = newVelZDR.transpose('Time','range','Vel')
#    data['Vel2ZeroZDR'] = newVelZDR
    return data

def globalAttr(data):
    data.attrs['Experiment']= 'TRIPEX-POL, Forschungszentrum Juelich'
    data.attrs['Instrument']= 'GRARAD-94 94.4 GHz cloud radar (W band)'
    data.attrs['Data']= 'Produced by Leonie von Terzi, lterzi@uni-koeln.de'
    data.attrs['Institution']= 'Data processed within the Emmy-Noether Group OPTIMIce, Institute for Geophysics and Meteorology, University of Cologne, Germany'
    data.attrs['Latitude']= '50.908547 N'
    data.attrs['Longitude']= '6.413536 E'
    data.attrs['Altitude']= 'Altitude of the JOYCE (www.joyce.cloud) platform: 111m asl'
    return data
def varAttr(data,spec=False,KDP=False):
    data.height.attrs['standard_name'] = 'height'
    data.height.attrs['long_name'] = 'projected height'
    data.height.attrs['units'] = 'm'

    if spec==True:
        data.HSpec.attrs['long_name'] = 'Doppler spectrum, horizontal polarization'
        data.HSpec.attrs['units'] = 'mm6 m-3'
        
        data.VSpec.attrs['long_name'] = 'Doppler spectrum, vertical polarization'
        data.VSpec.attrs['units'] = 'mm6 m-3'
        data.VSpec.attrs['comment'] = 'the variable ReVHSpec was saved incorrectly. Therefore, when the SNR is low, the spectral data at vertical resolution is not correct. In order to avoid using wrong data, you have to mask the VSpec, when the SNR_H is less than 10 dB!'     
 
        data.HSpec.attrs['long_name'] = 'Doppler spectrum, horizontal polarization'
        data.HSpec.attrs['units'] = 'mm6 m-3'

        data.HNoisePow.attrs['long_name'] = 'Integrated Doppler spectrum noise pover in horizontal channel'
        data.HNoisePow.attrs['units'] = 'mm6 m-3'

        data.VNoisePow.attrs['long_name'] = 'Integrated Doppler spectrum noise pover in vertical channel'
        data.VNoisePow.attrs['units'] = 'mm6 m-3'
        
        data.SNR_H.attrs['long_name'] = 'signal to noise ratio, horizontal polarization'
        data.SNR_H.attrs['units'] = 'dB'
        data.SNR_H.attrs['comment'] = 'the variable ReVHSpec was saved incorrectly. Therefore, when the SNR is low, the spectral data at vertical resolution is not correct. In order to avoid using wrong data, you have to mask the VSpec, when the SNR_H is less than 10 dB!'

        data.sSNR_H.attrs['long_name'] = 'signal to noise ratio spectrally resolved, horizontal polarization'
        data.sSNR_H.attrs['units'] = 'dB'
        data.sSNR_H.attrs['comment'] = 'the variable ReVHSpec was saved incorrectly. Therefore, when the SNR is low, the spectral data at vertical resolution is not correct. In order to avoid using wrong data, you have to mask the spectral bins of the VSpec, when the sSNR_H is less than 10 dB!'

        data.sSNR_V.attrs['long_name'] = 'signal to noise ratio spectrally resolved, vertical polarization'
        data.sSNR_V.attrs['units'] = 'dB'

        data.SNR_V.attrs['long_name'] = 'signal to noise ratio, vertical polarization'
        data.SNR_V.attrs['units'] = 'dB'

        data.ZH.attrs['standard_name'] = 'equivalent_reflectivity_factor'
        data.ZH.attrs['long_name'] = 'equivalent reflectivity factor, horizontal polarisation'
        data.ZH.attrs['units'] = 'mm6 m-3'

        data.ZV.attrs['standard_name'] = 'equivalent_reflectivity_factor'
        data.ZV.attrs['long_name'] = 'equivalent reflectivity factor, vertical polarisation'
        data.ZV.attrs['units'] = 'mm6 m-3'

        data.ZDR.attrs['standard_name'] = 'ZDR'
        data.ZDR.attrs['long_name'] = 'Differential reflectivity'
        data.ZDR.attrs['units'] = 'dB'

        data.sZDR.attrs['standard_name'] = 'spectral ZDR'
        data.sZDR.attrs['long_name'] = 'spectral Differential reflectivity'
        data.sZDR.attrs['units'] = 'dB'

        data.maxVelH.attrs['long_name'] = 'maximum Doppler velocity (spectral edge), horizontal polarization'
        data.maxVelH.attrs['units'] = 'ms-1'

        data.maxVelV.attrs['long_name'] = 'maximum Doppler velocity (spectral edge), vertical polarization'
        data.maxVelV.attrs['units'] = 'ms-1'

        data.maxVelZDR.attrs['long_name'] = 'maximum Doppler velocity (spectral edge), spectral ZDR'
        data.maxVelZDR.attrs['units'] = 'ms-1'

        data.Vel2ZeroH.attrs['comment'] = 'in lack of a better name, this variable is corrected for the offset in Doppler velocity, such that the slow spectral edge is set to 0 ms-1. Use this variable to plot the spectrum with the slowest edge at 0 ms-1. This is necessary, since the Doppler velocity for measurements at a non-nadir elevation is influenced by the horizontal wind velocity. Without this correction, the spectrum will look like a snake ;)'
        data.Vel2ZeroH.attrs['units'] = 'ms-1'  


    elif KDP==True:
        data.KDP_Alexander.attrs['standard_name'] = 'KDP'
        data.KDP_Alexander.attrs['long_name'] = 'Specific differential phase shift'
        data.KDP_Alexander.attrs['units'] = 'deg km-1'

    else:
        data.DBZ.attrs['standard_name'] = 'equivalent_reflectivity_factor'
        data.DBZ.attrs['units'] = 'dB'
 
        data.ZDR.attrs['standard_name'] = 'ZDR'
        data.ZDR.attrs['long_name'] = 'Differential reflectivity'
        data.ZDR.attrs['units'] = 'dB'

        data.KDP.attrs['standard_name'] = 'KDP'
        data.KDP.attrs['long_name'] = 'Specific differential phase shift'
        data.KDP.attrs['units'] = 'deg km-1'

        data.PhiDP.attrs['standard_name'] = 'PhiDP'
        data.PhiDP.attrs['long_name'] = 'Differential phase shift'
        data.PhiDP.attrs['units'] = 'deg'

        data.RHV.attrs['standard_name'] = 'rho_hv'
        data.RHV.attrs['long_name'] = 'Correlation coefficient'
        data.RHV.attrs['units'] = ''

        data.SK.attrs['standard_name'] = 'skewness'
        data.SK.attrs['long_name'] = 'skewness'
        data.SK.attrs['units'] = ''

        data.WIDTH.attrs['standard_name'] = 'doppler_spectrum_width'
        data.WIDTH.attrs['long_name'] = 'radar spectral width'
        data.WIDTH.attrs['units'] = 'ms-1'

        data.VEL.attrs['standard_name'] = 'radial_velocity_of_scatterers_away_from_instrument'
        data.VEL.attrs['long_name'] = 'mean Doppler velocity'
        data.VEL.attrs['units'] = 'ms-1'

        data.KURT.attrs['standard_name'] = 'kurtosis'
        data.KURT.attrs['long_name'] = 'kurtosis'
        data.KURT.attrs['units'] = ''
    return data
