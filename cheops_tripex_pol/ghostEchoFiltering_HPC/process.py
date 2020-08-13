###############################################################################
# now this script contains all the functions to deal with the W-band data, aswell as the saving to nc file.
# It is mainly from Jose, but I adjusted some things to deal better with xarray dataarrays and save time
###############################################################################
import os
import numpy as np
import pandas as pd
import xarray as xr
import glob as gb
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import scipy


def getNearestTimeIndex(refArr, target, tolerance):
# now I think this is not used, but need to check before removing..
    delta = np.abs(refArr - target)
    indexVal = np.arange(len(delta))
    indexVal = indexVal[delta <= tolerance]
    targetIndex = np.argmin(delta)
    
    try:
        testMin = np.min(indexVal)
        return testMin
        
    except:
        return None
    
    
def interpMaskToW(wData, kaData, maxVelKa, minVelKa):
#- now striclty speaking we dont need to interpolate the height from W-band to ka-band data, because of the regridding done, but for now I am leaving it because when I change something then for some reason it wont run anymore...
    rangeW = np.append(wData.C1Range.values,
                       wData.C2Range.values,)
    rangeW = np.append(rangeW, wData.C3Range.values,)
    rangeW = np.append(rangeW, wData.C4Range.values,)

    maxVelKaXR= xr.DataArray(maxVelKa, 
                             coords={'range':kaData.range},
                             dims={'range'})
    minVelKaXR= xr.DataArray(minVelKa, 
                             coords={'range':kaData.range}, 
                             dims={'range'})
    maxVelKaInt = maxVelKaXR.interp(range=rangeW)
    minVelKaInt = minVelKaXR.interp(range=rangeW)
    
    return maxVelKaInt, minVelKaInt


def spectraCleaning(wData, maxVelKa, minVelKa,channel):
# now mask the wband data. Jose had it all complicated, I didnt understand it so I just wrote my own too lines here ;) 
    maskedSpec = wData.where(wData.Vel > minVelKa)
    maskedSpec = maskedSpec.where(wData.Vel < maxVelKa)

    return maskedSpec

def writeOutPut(wTime, wIndex, outputPath, wSpec, kaNoise, kaSpec, channel, ZeKa, mdvKa, ZeW, mdvW,dv, ZeX=0, mdvX=0, xNoise=0, xSpec=0):
#- writes the spectra to ncfiles ( I am appending each timestep to an existing one, so If you want to generate new ncfiles, you need to delete the old ones. This is not ideal, so maybe later I can add a function which generated a new ncfile for the first second of each hour..)
# this is also mainly written by Jose, I just added some variables.

    kaNoise = kaNoise.reshape(1,len(kaSpec.range.values))
    #-- because we don't have co and cx channel in xband, we have to work around that, so only process if channel is co...
    if channel == 'o':
        xNoise = xNoise.reshape(1,len(xSpec.range.values))

    #### writing clean spectra ####
    idOutput = pd.to_datetime(wTime).strftime('%Y%m%d_%H')

    dataPath = os.path.join(outputPath, idOutput+'_LVL0_data_regridded_test.nc')

    try:
        rootgrp = Dataset(dataPath, 'a', format='NETCDF4')

    except:# IOError:
        rootgrp = creatNCFile(wSpec, dataPath, kaSpec, xSpec,dv)

    ###append new time
    rootgrp['time'][wIndex] = np.array(np.float(wTime)/10**9).reshape(1)
    if channel == 'x':
        idCh = 'H'
    if channel == 'o':
        idCh = 'V'
    #-- WData:
    rootgrp['WSpec'+idCh][wIndex] = wSpec['Spec'+idCh]
    rootgrp['W_Z_'+idCh][wIndex] = ZeW
    rootgrp['W_VEL_'+idCh][wIndex] = mdvW
    
    # also put noise level in nc file
    #rootgrp['WSpecNoise'+idCh][wIndex] = wSpec['Noise'+idCh]
     
    #-- KaData:
    rootgrp['KaSpec'+idCh][wIndex] = 10**(kaSpec/10)
    rootgrp['KaSpecNoise'+idCh][wIndex] = 10**(kaNoise/10)
    rootgrp['Ka_Z_'+idCh][wIndex] = ZeKa
    rootgrp['Ka_VEL_'+idCh][wIndex] = mdvKa


    #-- xData:
    if channel == 'o':
        rootgrp['XSpec'+idCh][wIndex] = 10**(xSpec/10)
        rootgrp['XSpecNoise'+idCh][wIndex] = 10**(xNoise/10)
        rootgrp['X_Z_'+idCh][wIndex] = ZeX
        rootgrp['X_VEL_'+idCh][wIndex] = mdvX

        
    rootgrp.close()


def creatNCFile(wSpec, dataPath, kaSpec, xSpec,dv):
    rootgrp = Dataset(dataPath, 'w', format='NETCDF4')
    rootgrp = globalAttributes(rootgrp)
    dv = dv.reshape(1)

    #dimension
    time = rootgrp.createDimension('time', None)

    Range = rootgrp.createDimension('range', kaSpec['range'].shape[0])
    #-- the velocity of Ka spectrum is called Doppler in the original nc file...
    DopplerKa = rootgrp.createDimension('dopplerKa', kaSpec['doppler'].shape[0])
    DopplerX = rootgrp.createDimension('dopplerX', xSpec['doppler'].shape[0])
    DopplerW =  rootgrp.createDimension('dopplerW', wSpec['Vel'].shape[0])
    dvW = rootgrp.createDimension('dvW', dv.shape[0])
    #dim variables 
    times = rootgrp.createVariable('time', np.int32,
                                   ('time',), zlib=True)
    times = timeAttributes(times)
    
    Range = rootgrp.createVariable('range', np.float32, 
                                      ('range'), zlib=True)
    Range = rangeAttributes(Range)
    DopplersKa = rootgrp.createVariable('dopplerKa', np.float32, 
                                    ('dopplerKa'), zlib=True)
    DopplersKa = velAttributes(DopplersKa) 
    DopplersX = rootgrp.createVariable('dopplerX', np.float32, 
                                    ('dopplerX'), zlib=True)
    DopplersX = velAttributes(DopplersX) 

    DopplersW = rootgrp.createVariable('dopplerW', np.float32, 
                                    ('dopplerW'), zlib=True)
    DopplersW = velAttributes(DopplersW) 

    #dv of the spectra of W
    dvW = rootgrp.createVariable('dvW', np.float32,
                                 ('dvW'), zlib=True)
    dvW = dvAttributes(dvW)
    #other variables (Spectra and Noise for H and V)
    WSpecV = rootgrp.createVariable('WSpecV', np.float32,
                                    ('time','range','dopplerW'),
                                    zlib=True)
    WSpecV = spectAttributesV(WSpecV)
    WSpecH = rootgrp.createVariable('WSpecH', np.float32,
                                    ('time','range','dopplerW'),
                                    zlib=True)
    WSpecH = spectAttributesH(WSpecH)
    # Noise
    WSpecNoiseV = rootgrp.createVariable('WSpecNoiseV', np.float32,
                                    ('time','range','dopplerW'),
                                    zlib=True)
    WSpecNoiseV = noiseAttributes(WSpecNoiseV)
    WSpecNoiseH = rootgrp.createVariable('WSpecNoiseH', np.float32,
                                    ('time','range','dopplerW'),
                                    zlib=True)
    WSpecNoiseH = noiseAttributes(WSpecNoiseH)
    
    #-- kaData:
    KaSpecV = rootgrp.createVariable('KaSpecV', np.float32,
                                    ('time','range','dopplerKa'),
                                    zlib=True)
    KaSpecV = spectAttributesV(KaSpecV)
    KaSpecH = rootgrp.createVariable('KaSpecH', np.float32,
                                    ('time','range','dopplerKa'),
                                    zlib=True)
    KaSpecH = spectAttributesH(KaSpecH)
    
    KaSpecNoiseH = rootgrp.createVariable('KaSpecNoiseH', np.float32,
                                    ('time','range'),
                                    zlib=True)
    KaSpecNoiseH = noiseAttributes(KaSpecNoiseH)
    KaSpecNoiseV = rootgrp.createVariable('KaSpecNoiseV', np.float32,
                                    ('time','range'),
                                    zlib=True)
    KaSpecNoiseV = noiseAttributes(KaSpecNoiseV)

    #-- xData:
    XSpecV = rootgrp.createVariable('XSpecV', np.float32,
                                    ('time','range','dopplerX'),
                                    zlib=True)
    XSpecV = spectAttributesV(XSpecV)
    #XSpecH = rootgrp.createVariable('XSpecH', np.float32,
    #                                ('time','range','dopplerX'),
    #                                zlib=True)
    #XSpecH = spectAttributesH(XSpecH)
    
    #XSpecNoiseH = rootgrp.createVariable('XSpecNoiseH', np.float32,
    #                                ('time','range'),
    #                                zlib=True)
    #XSpecNoiseH = noiseAttributes(XSpecNoiseH)
    XSpecNoiseV = rootgrp.createVariable('XSpecNoiseV', np.float32,
                                    ('time','range'),
                                    zlib=True)
    XSpecNoiseV = noiseAttributes(XSpecNoiseV)

    #-- Ze and mdv:
    XZeV = rootgrp.createVariable('X_Z_V', np.float32,
                                    ('time','range'),
                                    zlib=True)
    XZeV = ZeAttributesV(XZeV)
    
    KaZeV = rootgrp.createVariable('Ka_Z_V', np.float32,
                                    ('time','range'),
                                    zlib=True)
    KaZeV = ZeAttributesV(KaZeV)
    KaZeH = rootgrp.createVariable('Ka_Z_H', np.float32,
                                    ('time','range'),
                                    zlib=True)
    KaZeH = ZeAttributesH(KaZeH)

    WZeV = rootgrp.createVariable('W_Z_V', np.float32,
                                    ('time','range'),
                                    zlib=True)
    WZeV = VelAttributesV(WZeV)
    WZeH = rootgrp.createVariable('W_Z_H', np.float32,
                                    ('time','range'),
                                    zlib=True)
    WZeH = ZeAttributesH(WZeH)


    WmdvH = rootgrp.createVariable('W_VEL_H', np.float32,
                                    ('time','range'),
                                    zlib=True)
    WmdvH = VelAttributesH(WmdvH)
    WmdvV = rootgrp.createVariable('W_VEL_V', np.float32,
                                    ('time','range'),
                                    zlib=True)
    WmdvV = VelAttributesV(WmdvV)
    
    KamdvH = rootgrp.createVariable('Ka_VEL_H', np.float32,
                                    ('time','range'),
                                    zlib=True)
    KamdvH = VelAttributesH(KamdvH)
    KamdvV = rootgrp.createVariable('Ka_VEL_V', np.float32,
                                    ('time','range'),
                                    zlib=True)
    KamdvV = VelAttributesV(KamdvV)

    XmdvV = rootgrp.createVariable('X_VEL_V', np.float32,
                                    ('time','range'),
                                    zlib=True)
    XmdvV = VelAttributesV(XmdvV)


    #passing the values to dimension    
    DopplersW[:] = wSpec['Vel']
    
    Range[:] = kaSpec['range']
    DopplersKa[:] = kaSpec['doppler']
    
    DopplersX[:] = xSpec['doppler']
    dvW[:] = dv
    
    return rootgrp


def globalAttributes(rootgrpOut):
   
    rootgrpOut.Experiment = 'TRIPEX-POL, Forschungszentrum Juelich'
    rootgrpOut.Instrument = 'GRARAD-94 94.4 GHz cloud radar (W band)'
    rootgrpOut.chirp_number = '4'
    rootgrpOut.Data = 'Produced by Jose Dias, jdiasnet@uni-koeln.de and Leonie von Terzi, lterzi@uni-koeln.de'
    rootgrpOut.Institution = 'Data processed within the Emmy-Noether Group OPTIMIce, Institute for Geophysics and Meteorology, University of Cologne, Germany'
    rootgrpOut.Latitude = '50.908547 N'
    rootgrpOut.Longitude = '6.413536 E'
    rootgrpOut.Altitude = 'Altitude of the JOYCE (www.joyce.cloud) platform: 111m asl'
   
    return rootgrpOut


def timeAttributes(timeRef):
   
    timeRef.long_name = 'time in sec since 01.01.1970 00:00:00'
    timeRef.units = 's'

    return timeRef


def rangeAttributes(rangeRef):

    rangeRef.long_name = 'Range from antenna to the center of each range gate'
    rangeRef.units = 'm'

    return rangeRef


def velAttributes(vel):

    vel.long_name = 'Doppler velocity'
    vel.units = 'm s-1'

    return vel

def spectAttributesV(spec):
    
    spec.long_name = 'Doppler spectrum at vertical polarization'
    spec.units = 'mm6 m-3'
    return spec
def spectAttributesH(spec):

    spec.long_name = 'Doppler spectrum at horizontal polarization'
    spec.units = 'mm6 m-3'
    return spec

def noiseAttributes(noise):

    noise.long_name = 'spectral Noise level'
    noise.units = 'mm6 m-3' 
    return noise
def VelAttributesV(Vel):
    Vel.long_name = 'Mean Doppler Velocity, vertical polarization'
    Vel.units = 'm/s'
    return Vel
def VelAttributesH(Vel):
    Vel.long_name = 'Mean Doppler Velocity, horizontal polarization'
    Vel.units = 'm/s'
    return Vel
def ZeAttributesV(Ze):
    Ze.long_name = 'linear equivalent reflectivity factor, vertical polarization'
    Ze.units = 'm/s'
    return Ze
def ZeAttributesH(Ze):
    Ze.long_name = 'linear equivalent reflectivity factor, horizontal polarization'
    Ze.units = 'm/s'
    return Ze
def dvAttributes(dv):
    dv.long_name = 'doppler resolution, W'
    return dv
