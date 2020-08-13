import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import glob as gb
from sys import argv

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.optimize import curve_fit
data2test = 'Ka_DBZ_H'
def selTime(date,start,end,pathData):
    data= xr.open_dataset(pathData+date+'_tripex_pol_3fr_L1_mom.nc')
    Ze = data[[data2test,'W_DBZ_H']]
    #ZeW = data['W_DBZ_V']#.isel(dvW=0)
    return Ze.sel(time=slice(start,end))#,ZeKa.sel(time=slice(start,end))
def fun(x,B):
    return x + B

pathData = '/work/lvonterz/cheops_tripex_pol/output/'
Dates = ['20181201','20181215','20181223','20181227','20181227','20190110','20190122','20190131']
Starts = ['2018-12-01 05:00:00','2018-12-15 20:00:00','2018-12-23 06:00:00','2018-12-27 03:00:00','2018-12-27 16:00:00','2019-01-10 00:00:00','2019-01-22 00:00:00','2019-01-31 00:00:00']
Ends = ['2018-12-01 12:00:00','2018-12-15 23:59:59','2018-12-23 10:00:00','2018-12-27 14:00:00','2018-12-27 23:59:59','2019-01-10 14:00:00','2019-01-22 14:00:00','2019-01-31 10:00:00']
popt_all=[]
for date,start,end in zip(Dates,Starts,Ends):
    print(date);print(start);print(end)
    Ze = selTime(date,start,end,pathData)
    #if ind == 0:
    ZeAll = Ze
        #ZeKaall = ZeKa
    #else:
    #    ZeAll = xr.merge([ZeAll,Ze])
        #ZeKaall = xr.merge([ZeKaall,ZeKa])
    #ind += 1
    print(ZeAll)
    #only select Ze smaller -10
    ZeAll = ZeAll.where(ZeAll['W_DBZ_H']<=-10.0, np.NaN)
    ZeAll = ZeAll.where(ZeAll[data2test]<=-10.0, np.NaN)
    print(ZeAll[data2test].to_dataframe())
    ZeW = []
    ZeKa = []
    for tSel in ZeAll.time.values:
        ZeW.extend(ZeAll['W_DBZ_H'].sel(time=tSel).values)
        ZeKa.extend(ZeAll[data2test].sel(time=tSel).values)
        print(tSel)
    ZeW = np.array(ZeW)
    ZeKa = np.array(ZeKa)
    indW = ~np.isfinite(ZeW)
    indKa = ~np.isfinite(ZeKa)
    indBoth = indW + indKa
    ZeWFinite = ZeW[~indBoth]
    ZeKaFinite = ZeKa[~indBoth]
    popt,pcov = curve_fit(fun, ZeKaFinite,ZeWFinite)
#    plt.scatter(ZeKa,ZeW,c='b')
    plt.scatter(ZeKaFinite,ZeWFinite,c='r',alpha=0.5)
    plt.plot(ZeKaFinite,ZeKaFinite+popt,label='offset: %5.2f'%float(popt))
    plt.plot([-70,-5],[-70,-5])
    plt.xlim([-70,-5])
    plt.ylim([-70,-5])
    plt.xlabel(data2test)
    plt.ylabel('W_DBZ_H')
    plt.title(start+' '+end+'atm_att_corrected')
    plt.legend()
    plt.savefig(start+'_'+data2test+'.png')
    #plt.show()
    plt.close()
    print(popt)
    print(pcov)
    #quit()
    #plt.plot(
    popt_all.extend(popt)
print(popt_all)    
