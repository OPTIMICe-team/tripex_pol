import os
import numpy as np
import pandas as pd
import xarray as xr
import glob as gb
import matplotlib.pyplot as plt
from netCDF4 import Dataset



def loadWData(filePath):

    wData = xr.open_dataset(filePath)
    correctTime = wData.Time + wData.Timems*10**(-3)
    wData.Time.values = correctTime.values
    wData.Time.attrs['units']='seconds since 2001-01-01 00:00:00 UTC'
    wData = xr.decode_cf(wData)
    
    return wData
