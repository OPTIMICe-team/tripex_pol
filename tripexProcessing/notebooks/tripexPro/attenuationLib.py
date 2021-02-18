import atmFunc
import pyPamtra
import netCDF4
import numpy as np
import pandas as pd
import tripexLib as trLib
from scipy import interpolate


def changeAttListOrder(interpAttDataList, variable,
                        radarFreqs):
    varNames = variable.keys()
    attDFList = [0,0,0]
    
    for varName in varNames:
        radarIndex = varNames.index(varName)#dataFrameList
        freqIndex = radarFreqs.index(variable[varName]['freq'])#interpAttDataList
        attDFList[radarIndex] = interpAttDataList[freqIndex]
        
    return attDFList


def applyAttCorr(dataFrameList, interpAttDataList, variable):
    
    varNames = variable.keys()
    
    for varName in varNames:
    
        radarIndex = varNames.index(varName)#dataFrameList
        
        dataFrameList[radarIndex] = dataFrameList[radarIndex] + interpAttDataList[radarIndex]
        
    #print(dataFrameList)
    return dataFrameList


def convCloudTime2HumTime(time):
    hour = np.array(time, np.int)
    minuteDec = np.array((time - hour)*60)
    minute = np.array((time - hour)*60, np.int)
    secondDec = np.array((minuteDec - minute)*60)
    second = np.array((minuteDec - minute)*60, np.int)

    timeDF = pd.DataFrame()
    timeDF['hour'] = hour
    timeDF['minute'] = minute
    timeDF['second'] = np.round(secondDec, 6)
    timeDF['timeConv'] = timeDF.hour.astype(str)+':'+timeDF.minute.astype(str)+':'+timeDF.second.astype(str)
    
    return timeDF.timeConv


def getResampledTimeRange(rangeRef, rangeTolerance, timeRef,
                          time, timeTolerance, dataToResample,
                          year, month, day, height_M):
   
    year = int(year)
    month = int(month)
    day = int(day)    

    usedIndexRange = np.ones(len(rangeRef))*np.nan
    timeConv = convCloudTime2HumTime(time)
    humanTime = pd.datetime(year, month, day) + pd.to_timedelta(np.array(timeConv), unit='s')
    #humanTime=time
    #Remove duplicates
    varData = pd.DataFrame()
    varDataTemp = pd.DataFrame(index=humanTime, columns=height_M, data=dataToResample)
    varDataTemp['times'] = humanTime
    varDataTemp = varDataTemp.sort_values(by=['times'], ascending=[True])
    varDataTemp = varDataTemp.drop_duplicates(subset=['times'])
    varData = varData.append(varDataTemp)

    #Resample in time
    timeIndexList = trLib.getIndexList(varData, timeRef, timeTolerance)
    emptyDataFrame = pd.DataFrame(index=timeRef, columns=varData.columns)
    resampledTime = trLib.getResampledDataPd(emptyDataFrame, varData, timeIndexList)
    timeDeviation = trLib.getDeviationPd(timeRef, resampledTime, timeTolerance)
    del resampledTime['times']
    resampledTime = resampledTime.T

    #Resample in range
    resampledTime['ranges'] = height_M
    rangeIndexList = trLib.getIndexList(resampledTime, rangeRef, rangeTolerance)
    emptyDataFrame = pd.DataFrame(index=rangeRef, columns=resampledTime.columns)
    resampledTimeRange = trLib.getResampledDataPd(emptyDataFrame, resampledTime,
                                             rangeIndexList)
    rangeDeviation = trLib.getDeviationPd(rangeRef, resampledTimeRange)
    del resampledTimeRange['ranges']

    resampledTimeRangeArr = np.array(resampledTimeRange[:],float)
    resampledTimeRangeArr = np.ma.masked_invalid(resampledTimeRangeArr)

    return resampledTimeRangeArr


def getInterpData(time, timeRef, height_M,
                  resampledData, 
                  dataToInterp,
                  rangeRef):
    
    funcToInterp2d = interpolate.interp2d(time, height_M, dataToInterp.T)

    secTimeRef = timeRef.hour + timeRef.minute/(60.) + timeRef.second/(60.*60.)
    interpData = funcToInterp2d(secTimeRef, rangeRef)
    
    qualityFlag = np.zeros_like(resampledData, np.float)
    interpData[~resampledData.mask] = resampledData[~resampledData.mask]
    qualityFlag[~resampledData.mask] = 1

    return interpData, qualityFlag


def getDescriptor():
    
    descriptorFile = np.array([
      #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as'
      # 'beta_as' # 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 
      #'d_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
      ('cwc_q', -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0, 3, 1, 
        'mono', -99.0, -99.0, -99.0, -99.0, 2e-05, -99.0, 'mie-sphere',
        'khvorostyanov01_drops', -99.0)], 
      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), 
             ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), 
             ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), 
             ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), 
             ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), 
             ('d_2', '<f8'), ('scat_name', 'S15'), ('vel_size_mod', 'S30'), 
             ('canting', '<f8')])
    
    return descriptorFile 


def getAtmAttPantra(cloudNetFilePath, radarFreqs):

    cloudNetData = netCDF4.Dataset(cloudNetFilePath,'r')
    height_M = cloudNetData['model_height'][:]-105# added -105[m] here to Joses original file since the cloudnet mocel_height start at msl and radar data is at 105 above msl
    temp = cloudNetData['temperature'][:]#[k]
    speHum = cloudNetData['specific_humidity'][:]#[kg/kg]
    press = cloudNetData['pressure'][:]#[Pa]
    time = cloudNetData['time'][:]#[hour fraction]
    att_atm_cloud = cloudNetData['radar_gas_atten'][:]#[dB]
    att_liq_cloud = cloudNetData['radar_liquid_atten'][:]#[dB]

    relHum = atmFunc.speHumi2RelHum(speHum, temp, press)#[%]
    vaporPress = atmFunc.calcVaporPress(speHum, temp, press)#[Pa] (!! variable name)
    waterVaporDens = atmFunc.calcVaporDens(vaporPress, temp)#[kg/m^3]
    dryAirDens = atmFunc.calcDryAirDens(press, waterVaporDens, temp)#[kg/m^3]
    
    descriptorFile = getDescriptor()

    pam = pyPamtra.pyPamtra()
    for hyd in descriptorFile:
        pam.df.addHydrometeor(hyd)

    pam.nmlSet['active'] = True
    pam.nmlSet['passive'] = False
    pam.nmlSet["radar_attenuation"] = 'bottom-up'

    pamData = dict()
    heightArr= np.zeros((len(time), len(height_M))) + height_M
    pamData['hgt'] = heightArr
    pamData['temp'] = temp
    pamData['relhum'] = relHum
    pamData['press'] = press

    pam.createProfile(**pamData)
    pam.runParallelPamtra(radarFreqs, pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=8)

    return pam.r, time, height_M, temp, relHum, press


def getInterpQualFlagList(results, time, timeRef, timeTolerance,
                          height_M, rangeRef, rangeTolerance,
                          radarFreqs, year, month, day):
    gridy = 0
    interpDataList = []
    qualityFlagList = []

    for i, radarFreq in enumerate(radarFreqs):
        interpDataList.append(pd.DataFrame())
        qualityFlagList.append(pd.DataFrame())

    for i, radarFreq in enumerate(radarFreqs):
    
        attAtm=results['Att_atmo'][:,gridy,:, i]
        #attHydro=results['Att_hydro'][:,gridy,:, i]
    
        attAtm2Way = 2*np.cumsum(attAtm,axis=1)
#        print(attAtm2Way)
        resampledData = getResampledTimeRange(rangeRef, rangeTolerance, timeRef,
                                              time, timeTolerance, attAtm2Way,
                                              year, month, day, height_M)

        interpData, qualityFlag = getInterpData(time, timeRef, height_M,
                                                resampledData, attAtm2Way,
                                                rangeRef)
    
        interpDF = pd.DataFrame(index=timeRef, columns=rangeRef,
                                data=interpData.T)
    
        qualityFlagDF = pd.DataFrame(index=timeRef, columns=rangeRef,
                                     data=qualityFlag.T)

        interpDataList[i] = interpDataList[i].append(interpDF)
        qualityFlagList[i] = qualityFlagList[i].append(qualityFlagDF)
        
    return interpDataList, qualityFlagList
