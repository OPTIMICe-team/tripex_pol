import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import xarray as xr
from sys import argv
import glob
import os
from netCDF4 import Dataset
import time
def regridKaData(date, files,variableName,timeFreq,rangeRef,
                timeTol,rangeTol,rangeOffset):
    print(date)
    Data = xr.open_dataset(files)
    Data.time.values = Data.time.values + Data.microsec.values*10**(-6) #the file stores time a bit weird... also the units need to be changed, in order for xarray to have the correct time axis
    Data.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
    Data = xr.decode_cf(Data)
    #-- in order to force the dataset to start at 00:00:00 (resample doesnt have the option, and reindex needs a lot of memory...) I check whether the first time is already 00:00:00, and if not I force it to, and will later remove the values for the first time in the dataset
    Datastart = pd.to_datetime(Data.time.values[0])
    DataStartStr = Datastart.strftime('%S')
    DataHour = Datastart.strftime('%H')
    DataStart = '00'
    Dataend = pd.to_datetime(Data.time.values[-1])
    DataEndStr = Dataend.strftime('%M')
    DataEnd = '59'
    print('DataStart:',DataStart)
    print('DataHour:',DataHour)
    print('DataStartData:',Datastart)
    
    if DataStart != DataStartStr:
        Data.time.values[0] = pd.to_datetime(date + pd.offsets.Hour(int(DataHour)))
        mask = True
    else:
        mask = False
    print('mask:',mask)
    print(Data.time)
    if DataEnd != DataEndStr:
        Data.time.values[-1] = pd.to_datetime(date) + pd.offsets.Hour(int(DataHour)) + pd.offsets.Hour(1) + pd.offsets.Minute(59) + pd.offsets.Second(59)
        maskLast = True
    else:
        maskLast = False
    print(Data.time)
    print('mask:',mask)
    #-- since the first range gate of the Ka-band is not in the same height as the one from w band, we have to correct for that
    Data.range.values = Data.range.values + rangeOffset

    #-- resample the data to the reference time and height grid, using the tolerance defined before 
    #-- get rid of duplicate time indices...
    _, index_time = np.unique(Data['time'], return_index=True)
    Data = Data.isel(time=index_time)
    DataRegridded = Data[['RadarConst','elv','npw1']].resample(time=timeFreq).nearest(tolerance=timeTol)
    for ind,variable in enumerate(variableName):
        var = [variable]
        print(var)
        print(Data[var])
        #-- resample time to the timefrequency we want (it's faster than reindex, and we have a lot of timestepts...))
        DatanewTime = Data[var].resample(time=timeFreq).nearest(tolerance=timeTol)
        #-- resample range (in theory you could resample time and range together with reindex, but that takes a lot of memory)
        DatanewRange = DatanewTime.reindex({'range':rangeRef}, method='nearest',tolerance=rangeTol)
        DataRegriddedNew = DatanewRange.where(DataRegridded.elv == 90)
        if mask == True:
            DataRegriddedNew = DataRegriddedNew.where(DataRegriddedNew.time != Data.time.values[0], np.NaN)
        if maskLast == True:
                DataRegriddedNew = DataRegriddedNew.where(DataRegriddedNew.time != Data.time.values[-1], np.NaN)
        DataRegridded = xr.merge([DataRegridded,DataRegriddedNew])
    return DataRegridded, DataHour


def regridWData(date, files,variableName,timeFreq,
                rangeRefC1,rangeRefC2,rangeRefC3,rangeRefC4,
                timeTol,rangeTol):
    Data = xr.open_dataset(files)
    Data.Time.values = Data.Time.values + Data.Timems*10**(-3)
    Data.Time.attrs['units']='seconds since 2001-01-01 00:00:00 UTC'
    Data = xr.decode_cf(Data)
    ChirpNum = Data.ChirpNum.values
    print(Data)
    #-- in order to force the dataset to start at 00:00:00 (resample doesnt have the option, and reindex needs a lot of memory...) I check whether the first time is already 00:00:00, and if not I force it to, and will later remove the values for the first time in the dataset
    Datastart = pd.to_datetime(Data.Time.values[0])
    DataStartStr = Datastart.strftime('%S')
    DataHour = Datastart.strftime('%H')
    DataStart = '00'
    Dataend = pd.to_datetime(Data.Time.values[-1])
    DataEndStr = Dataend.strftime('%M')
    DataEnd = '59'
    print('DataStart:',DataStartStr)
    if DataStart != DataStartStr:
        Data.Time.values[0] = pd.to_datetime(date + pd.offsets.Hour(int(DataHour)))
        mask = True
    else:
        mask = False
    print('mask:',mask)
    if DataEnd != DataEndStr:
        Data.Time.values[-1] = pd.to_datetime(date + pd.offsets.Hour(int(DataHour)) + pd.offsets.Minute(59) + pd.offsets.Second(59))
        maskLast = True
    else:
        maskLast = False
    print(Data)
    print('mask:',mask)
    #-- we need to transfer dataset into pandas dataframe in order to get rid of duplicate time indices...Dataframe = Data[variableName].to_dataframe()
    _, index_time = np.unique(Data['Time'], return_index=True)
    Data = Data.isel(Time=index_time)
    DataRegridded = Data[variableName].resample(Time=timeFreq).nearest(tolerance=timeTol)
    #-- because the dopplerlen, maxvel and elv are not dependent on range, we only need to resample them into the correct time...
    for variable in ['HSpec','VSpec']:#,'HNoisePow','VNoisePow']:
        print(variable)
        for ic in range(ChirpNum):        
            chirp = 'C{index}'.format(index=ic+1)
            variableName1=[chirp+variable]
            print('Databefore:',Data[variableName1].values)
            if ic == 0:
                rangeRef = rangeRefC1
            elif ic == 1:
                rangeRef = rangeRefC2
            elif ic == 2:
                rangeRef = rangeRefC3
            else:
                rangeRef = rangeRefC4

            ranges = chirp+'Range'
            print(ranges)
            #-- resample the data to the reference time and height grid, using the tolerance defined before
           #            DatanewTime = Dataset.reindex({'Time':timeRef},method='nearest',tolerance=timeTol)
            DatanewTime = Data[variableName1].resample(Time=timeFreq).nearest(tolerance=timeTol)
            DataRegriddedNew = DatanewTime.reindex({ranges:rangeRef}, method='nearest',tolerance=rangeTol)
            if mask == True:
                DataRegriddedNew = DataRegriddedNew.where(DataRegriddedNew.Time != Data.Time.values[0], np.NaN)
            if maskLast == True:
                DataRegriddedNew = DataRegriddedNew.where(DataRegriddedNew.Time != Data.Time.values[-1], np.NaN)
            DataRegridded = xr.merge([DataRegridded,DataRegriddedNew])
            print(DataRegridded[variableName1].values)
    return DataRegridded, DataHour


def regridXData(date, files,variableName,timeFreq,rangeRef,
                timeTol,rangeTol,rangeOffset):
    Data = xr.open_dataset(files)
    Data.time.values = Data.time.values + Data.microsec.values*10**(-6) #the file stores time a bit weird... also the units need to be changed, in order for xarray to have the correct time axis
    Data.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
    Data = xr.decode_cf(Data)
    #-- in order to force the dataset to start at 00:00:00 (resample doesnt have the option, and reindex needs a lot of memory...) I check whether the first time is already 00:00:00, and if not I force it to, and will later remove the values for the first time in the dataset
    
    Datastart = pd.to_datetime(Data.time.values[0])
    DataStartStr = Datastart.strftime('%S')
    DataHour = Datastart.strftime('%H')
    DataDay = Datastart.strftime('%Y%m%d %H:00:00')
    DataEnd = Datastart+pd.offsets.Hour(1)
    DataDayEnd = DataEnd.strftime('%Y%m%d %H:59:56')
    DataStart = '00'
    timeRef = pd.date_range(DataDay,DataDayEnd, freq = timeFreq)
    #-- since the first range gate of the X-band is not in the same height as the one from w band, we have to correct for that
    Data.range.values = Data.range.values + rangeOffset
    '''
    print('DataStart:',DataStartStr)
    if DataStart != DataStartStr:
        Data.time.values[0] = pd.to_datetime(date + pd.offsets.Hour(int(DataHour)))
        mask = True
    else:
        mask = False
    print('mask:',mask)
    '''
    #-- resample the data to the reference time and height grid, using the tolerance defined before 
    #-- get rid of duplicate time indices...
    _, index_time = np.unique(Data['time'], return_index=True)
    Data = Data.isel(time=index_time)
    DataRegridded = Data[['RadarConst','elv','npw1']].reindex({'time':timeRef}, method='nearest',tolerance=timeTol)
    for ind,var in enumerate(variableName):
        print(var)
        print(Data[var])
        #-- since the x-band has so many doppler values, I have to resample differently, otherwise i am running out of memory
        if (var == 'SPCco') or (var == 'SPCcx'):
            DataDoppler = Data.doppler.where(Data.doppler >= -20.0, drop=True)
            DataDoppler = DataDoppler.where(DataDoppler <= 20.0, drop=True)
            dataEmpty = np.empty((len(timeRef),len(rangeRef),len(DataDoppler.values)))
            Dataregriddednew = xr.DataArray(dataEmpty, dims=('time','range','doppler'),
                          coords={'time':timeRef,'range':rangeRef, 'doppler':DataDoppler})
            print(timeRef)
            print(timeTol)
            for tind,ti in enumerate(timeRef.values):
                print(ti)
                try:
                    DataStep = Data[var].sel(time=ti, method='nearest',tolerance=timeTol)
                except:
                    DataStep = Dataregriddednew[tind,:,:]
                DataStep = DataStep.where(DataStep.doppler >= -20.0, drop=True)
                DataStep = DataStep.where(DataStep.doppler <= 20.0, drop=True)
                DataStep = DataStep.reindex({'range':rangeRef}, method='nearest',tolerance=rangeTol)
                Dataregriddednew[tind,:,:] = DataStep
            Dataregriddednew.name = var
            
        else:
            DatanewDoppler = Data[var]
            Dataregriddednew = DatanewDoppler.reindex({'range':rangeRef}, method='nearest',tolerance=rangeTol)
            Dataregriddednew = Dataregriddednew.reindex({'time':timeRef}, method='nearest',tolerance=timeTol)
        print('Dataregriddednew',Dataregriddednew)
        DataRegridded = xr.merge([DataRegridded,Dataregriddednew])
    return DataRegridded, DataHour

























