import glob
import numpy as np
import pandas as pd

import xarray as xr
from sys import argv
import glob
import os
from netCDF4 import Dataset
import time
'''
def regridKaXData(files,variableName,timeFreq,rangeRef,
                timeTol,rangeTol):
    print(variableName)
    Data = xr.open_dataset(files)
    Data.time.values = Data.time.values + Data.microsec.values*10**(-6) #the file stores time a bit weird... also the units need to be changed, in order for xarray to have the correct time axis
    Data.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
    Data = xr.decode_cf(Data)
    #-- resample the data to the reference time and height grid, using the tolerance defined before 
    #--  get rid of duplicate time indices...
    index_time = np.unique(Data['time'], return_index=True)
    Data = Data.isel(time=index_time)

    for ind,variable in enumerate(variableName):
        if ind == 0:
            var = [variable,'elv']
        else:
            var = [variable,'elv']
        print(var)
        DatanewTime = Data[var].resample(time=timeFreq).nearest(tolerance=timeTol)
        print(DatanewTime)
        if ind == 0:
            DatanewRange = DatanewTime.reindex({'range':rangeRef}, method='nearest', tolerance=rangeTol)
            DataRegridded = DatanewRange.where(DatanewRange.elv == 90)
        else:
            DatanewRange = DatanewTime.reindex({'range':rangeRef}, method='nearest',tolerance=rangeTol)
            DataRegriddedNew = DatanewRange.where(DatanewRange.elv == 90)
            DataRegridded = xr.merge([DataRegridded,DataRegriddedNew])
    return DataRegridded
'''
def regridKaXData(files,variableName,timeFreq,rangeRef,
                timeTolerance,rangeTolerance):
    print(variableName)
    Data = xr.open_dataset(files)
    Data.time.values = Data.time.values + Data.microsec.values*10**(-6) #the file stores time a bit weird... also the units need to be changed, in order for xarray to have the correct time axis
    Data.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
    Data = xr.decode_cf(Data)
    #-- resample the data to the reference time and height grid, using the tolerance defined before 
    #-- get rid of duplicate time indices...
    _, index_time = np.unique(Data['time'], return_index=True)
    Data = Data.isel(time=index_time)
    for ind,variable in enumerate(variableName):
        if ind == 0:
            var = [variable,'elv']
        else:
            var = [variable,'elv']
        print(var)
        DatanewTime = Data[var].resample(time=timeFreq).nearest(tolerance=timeTolerance)
        print(DatanewTime)
        if ind == 0:
            DatanewRange = DatanewTime.reindex({'range':rangeRef}, method='nearest')
            DataRegridded = DatanewRange.where(DatanewRange.elv == 90)
        else:
            DatanewRange = DatanewTime.reindex({'range':rangeRef}, method='nearest')
            DataRegriddedNew = DatanewRange.where(DatanewRange.elv == 90)
            DataRegridded = xr.merge([DataRegridded,DataRegriddedNew])
    return DataRegridded

#-- regrid the Wband data, this is a bit more complicated because we have different chirps with different range frequencies...
def regridWData(files,variableName,timeRef,rangeRef,
                timeTolerance,rangeTolerance):
    Data = xr.open_dataset(files)
    Data.Time.values = Data.Time.values + Data.Timems*10**(-3)
    Data.Time.attrs['units']='seconds since 2001-01-01 00:00:00 UTC'
    Data = xr.decode_cf(Data)
    #-- we need to transfer dataset into pandas dataframe in order to get rid of duplicate time indices...Dataframe = Data[variableName].to_dataframe()
    Dataframe = Data[variableName].to_dataframe()
    Dataframe = Dataframe.loc[~Dataframe.index.duplicated(keep='first')]
    Dataset = Dataframe.to_xarray()
    print(Dataset)
    DataRegridded = Dataset.reindex({'time':timeRef},method='nearest',tolerance=timeTol)
    #-- because the dopplerlen, maxvel and elv are not dependent on range, we only need to resample them into the correct time...
    for variable in ['VNoisePow','HNoisePow','VSpec','HSpec']:
        xr.concat([Data['C1'+variable], Data['C2'+variable], Data['C3'+variable],Data['C4'+variable]])
        chirp = 'C{index}'.format(index=ic+1)
        variableName1=[chirp+variable]
        #-- resample the data to the reference time and height grid, using the tolerance defined before
        Dataframe = Data[variableName1].to_dataframe()
        Dataframe = Dataframe.loc[~Dataframe.index.duplicated(keep='first')]
        Dataset = Dataframe.to_xarray()
        DatanewTime = Dataset.resample(Time=timeFreq).nearest(tolerance=timeTolerance)
        DataRegriddedNew = DatanewTime.reindex({ranges:rangeRef}, method='nearest')
        DataRegridded = xr.merge([DataRegridded,DataRegriddedNew])
    return DataRegridded



'''
import glob
import numpy as np
import pandas as pd

import xarray as xr
from sys import argv
import glob
import os
from netCDF4 import Dataset
import time

def regridKaXData(files,variableName,timeFreq,rangeRef,
                timeTolerance,rangeTolerance):
    print(variableName)
    Data = xr.open_dataset(files)
    Data.time.values = Data.time.values + Data.microsec.values*10**(-6) #the file stores time a bit weird... also the units need to be changed, in order for xarray to have the correct time axis
    Data.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
    Data = xr.decode_cf(Data)
    #-- resample the data to the reference time and height grid, using the tolerance defined before 
    #-- we need to transfer dataset into pandas dataframe in order to get rid of duplicate time indices...
    for ind,variable in enumerate(variableName):
        if ind == 0:
            var = [variable,'elv']
        else:
            var = [variable,'elv']
        print(var)
        DatanewRange = Data.reindex({'range':rangeRef}, method='nearest')
        Dataframe = DatanewRange[var].to_dataframe()
        Dataframe = Dataframe.loc[~Dataframe.index.duplicated(keep='first')]
        Dataset = Dataframe.to_xarray()
        if ind == 0:
            DataRegridded = Dataset.resample(time=timeFreq).nearest(tolerance=timeTolerance)
        #if ind == 0:
        #    DataRegridded = DatanewTime.reindex({'range':rangeRef}, method='nearest')
            #-- now we also have to mask data where the elevation at ka-band wasn't 90°:
            DataRegridded = DataRegridded.where(DataRegridded.elv == 90)
        else:
            DataRegriddedNew = Dataset.resample(time=timeFreq).nearest(tolerance=timeTolerance)
            DataRegriddedNew = DataRegriddedNew.where(DataRegriddedNew.elv == 90)
            DataRegridded = xr.merge([DataRegridded,DataRegriddedNew])
    return DataRegridded

#-- regrid the Wband data, this is a bit more complicated because we have different chirps with different range frequencies...
def regridWData(files,variableName,timeRef,rangeRef,
                timeTolerance,rangeTolerance):
    Data = xr.open_dataset(files)
    Data.Time.values = Data.Time.values + Data.Timems*10**(-3)
    Data.Time.attrs['units']='seconds since 2001-01-01 00:00:00 UTC'
    Data = xr.decode_cf(Data)
    #-- we need to transfer dataset into pandas dataframe in order to get rid of duplicate time indices...Dataframe = Data[variableName].to_dataframe()
    Dataframe = Data[variableName].to_dataframe()
    Dataframe = Dataframe.loc[~Dataframe.index.duplicated(keep='first')]
    Dataset = Dataframe.to_xarray()
    print(Dataset)
    DataRegridded = Dataset.reindex({'time':timeRef},method='nearest',tolerance=timeTol)
    #-- because the dopplerlen, maxvel and elv are not dependent on range, we only need to resample them into the correct time...
    for variable in ['VNoisePow','HNoisePow','VSpec','HSpec']:
        xr.concat([Data['C1'+variable], Data['C2'+variable], Data['C3'+variable],Data['C4'+variable]])
        chirp = 'C{index}'.format(index=ic+1)
        variableName1=[chirp+variable]
        #-- resample the data to the reference time and height grid, using the tolerance defined before
        Dataframe = Data[variableName1].to_dataframe()
        Dataframe = Dataframe.loc[~Dataframe.index.duplicated(keep='first')]
        Dataset = Dataframe.to_xarray()
        DatanewTime = Dataset.resample(Time=timeFreq).nearest(tolerance=timeTolerance)
        DataRegriddedNew = DatanewTime.reindex({ranges:rangeRef}, method='nearest')
        DataRegridded = xr.merge([DataRegridded,DataRegriddedNew])
    return DataRegridded

'''




























