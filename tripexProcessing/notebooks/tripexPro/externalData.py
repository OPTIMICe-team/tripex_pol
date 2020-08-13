import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import tripexLib as trLib
import glob

def getDFResampledTime(timeRef, timeTolerance, 
                       dataFrameToResample):
    
    #Remove duplicates
    varDataTemp = dataFrameToResample
    varDataTemp['times'] = varDataTemp.index
    varDataTemp = varDataTemp.sort_values(by=['times'], ascending=[True])
    varDataTemp = varDataTemp.drop_duplicates(subset=['times'])
    

    #Resample in time
    timeIndexList = trLib.getIndexList(varDataTemp, timeRef, timeTolerance)
    emptyDataFrame = pd.DataFrame(index=timeRef, columns=varDataTemp.columns)
    resampledTime = trLib.getResampledDataPd(emptyDataFrame, varDataTemp, timeIndexList)
    timeDeviation = trLib.getDeviationPd(timeRef, resampledTime, timeTolerance)
    del resampledTime['times']
    #resampledTime = resampledTime.T

    return resampledTime


def getDataRadiometer(year, month, day, timeRef,
                      timeTolerance, variable):

    if variable == 'LWP':
        fileID = 'sups_joy_mwr00_l2_clwvi_v01_'
        varName = 'clwvi'
    if variable == 'IWV':
        fileID = 'sups_joy_mwr00_l2_prw_v01_'
        varName = 'prw'
        
    radioPath = '/data/hatpro/jue/data/level2'    
    radioFile = ('/').join([radioPath, year[-2:]+month,
                            fileID+year+month+day+'*nc'])
    
    radioFileList=glob.glob(radioFile)
    
    radioFile = Dataset(radioFileList[0])
    radioTime = radioFile.variables['time'][:]
    varData = radioFile.variables[varName][:]
    flag = radioFile.variables['flag'][:]

    epoch = pd.datetime(1970, 1, 1)
    times = epoch + pd.to_timedelta(radioTime, 's')
    
    variableDF = pd.DataFrame(index=times, data=varData,
                              columns=[variable])
    variableDF['flag'] = flag
    
    resampledVariableDF = getDFResampledTime(timeRef, timeTolerance,
                                             variableDF)
    
    return resampledVariableDF


def getDataPluvio(year, month, day, timeRef,
                  timeTolerance):
    
    pluvioPath = '/data/hatpro/jue/data/pluvio/netcdf'
    pluvioIDFile = 'pluvio2_jue_'

    pluvioFile = ('/').join([pluvioPath, year[-2:]+month,
                             pluvioIDFile+year+month+day+'.nc'])

    pluvioData = Dataset(pluvioFile)
    timePl = pluvioData.variables['time'][:]
    rain_RT = pluvioData.variables['rain_rate'][:]
    accum_RT = pluvioData.variables['r_accum_RT'][:]
    accum_NRT = pluvioData.variables['r_accum_NRT'][:]
    tot_accum_NRT = pluvioData.variables['total_accum_NRT'][:]
    tot_accum_NRT = tot_accum_NRT-tot_accum_NRT[0]
    
    epoch = pd.datetime(1970, 1, 1)
    timesPl = epoch + pd.to_timedelta(timePl, 's')
    
    pluvioDF = pd.DataFrame(index=timesPl)
    pluvioDF['times'] = timePl
    pluvioDF['rainRT'] = rain_RT
    pluvioDF['accmRT'] = accum_RT
    pluvioDF['accumNRT'] = accum_NRT*60.
    pluvioDF['totAccumNRT'] = tot_accum_NRT
    
    pluvioDF = pluvioDF[pluvioDF.times > 0]
    del pluvioDF['times']
    
    resampledPluvioDF = getDFResampledTime(timeRef, timeTolerance,
                                           pluvioDF)
    
    return resampledPluvioDF


def getDataCeilo(year, month, day, timeRef,
                 timeTolerance):
    
    ceiloPath = '/data/TR32/D2/data/ceilo/level0b'
    ceiloIDFile = '_ct25k_jue_l0b.nc'

    ceiloFilePath = ('/').join([ceiloPath,year[-2:]+month,
                               year[-2:]+month+day+ceiloIDFile])
    
    ceiloData = Dataset(ceiloFilePath)
    ceiloBaseTime = ceiloData.variables['base_time'][:]
    ceiloTime = ceiloData.variables['time'][:]
    cloudBaseHeight = ceiloData.variables['first_cbh'][:]
    
    ceiloTime = ceiloTime*60*60 + ceiloBaseTime
    epoch = pd.datetime(1970, 1, 1)
    ceiloTime = epoch + pd.to_timedelta(ceiloTime, 's')
    
    ceiloDF = pd.DataFrame(index=ceiloTime,
                       data=cloudBaseHeight,
                       columns=['cloudBaseHeight'])
    
    ceiloDF = ceiloDF[ceiloDF.cloudBaseHeight > 0]
    
    resampledCeiloDF = getDFResampledTime(timeRef, timeTolerance,
                                          ceiloDF)
    
    return resampledCeiloDF
