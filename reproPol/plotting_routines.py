# this skript is a collection of the plotting routines required for plotting the tripex-pol dataset correctly
# Author: Leonie von Terzi


import os
import time
import xarray as xr
import numpy as np
import netCDF4 as nc
import glob as glob
import pandas as pd
import matplotlib.pyplot as plt

def plotMom(dataset, plotOutPath, strDate, fileID, dbz=False, vel=False,width=False,sk=False):
# plot the non-polarimetric moments
    fig, axes = plt.subplots(nrows=3, figsize=(18,18))
    if width==True:
        varX = 'X_WIDTH_H'
        varW = 'W_WIDTH_H'
        varKa = 'Ka_WIDTH_H'
        lim = (0,1)
    elif sk==True:
        varX = 'X_SK_H'
        varW = 'W_SK_H'
        varKa = 'Ka_SK_H'
        lim = (-1,1)
    elif dbz==True:
        varX = 'X_DBZ_H'
        varW = 'W_DBZ_H'
        varKa = 'Ka_DBZ_H'
        lim = (-35,25)
    elif vel==True:
        varX = 'X_VEL_H'
        varW = 'W_VEL_H'
        varKa = 'Ka_VEL_H'
        lim = (-3,1)

    radData = {varX:{'data':dataset, 'axis':axes[0], 'lim':lim},
               varKa:{'data':dataset, 'axis':axes[1], 'lim':lim},
               varW:{'data':dataset, 'axis':axes[2], 'lim':lim},
                }

    for rad in radData.keys():

        plot = radData[rad]['data'][rad].T.plot(ax=radData[rad]['axis'],
                                           vmax=radData[rad]['lim'][1],
                                           vmin=radData[rad]['lim'][0],
                                           cmap='jet')

        radData[rad]['axis'].set_title(rad +' '+ strDate)
        plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
        radData[rad]['axis'].grid()
        radData[rad]['axis'].set_xlabel('')

    plotFileName = ('_').join([strDate,fileID])
    filePathName = ('/').join([plotOutPath,plotFileName])
    plt.savefig(filePathName+'.png',dpi=200,bbox_inches='tight')

def plotDWR(DWR_X_Ka, DWR_Ka_W, plotOutPath, strDate, plotID):
    #print(dataset)
    fig, axes = plt.subplots(nrows=2, figsize=(18,12))
    radData = {'DWR_X_Ka':{'data':DWR_X_Ka, 'axis':axes[0], 'lim':(-5,20)},
               'DWR_Ka_W':{'data':DWR_Ka_W, 'axis':axes[1], 'lim':(-5,20)}}
    for rad in radData.keys():

        plot = radData[rad]['data'].T.plot(ax=radData[rad]['axis'],
                                        vmax=radData[rad]['lim'][1],
                                        vmin=radData[rad]['lim'][0],
                                        cmap='jet')

        radData[rad]['axis'].set_title(rad +' '+ strDate)
        plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
        radData[rad]['axis'].grid()
        radData[rad]['axis'].set_xlabel('')

    plotFileName = ('_').join([strDate,plotID])
    filePathName = ('/').join([plotOutPath,plotFileName])
    plt.savefig(filePathName+'.png',dpi=200, bbox_inches='tight')

def plotSpectra(dataset,dataLVL2,timeseries, plotOutPath, plotId, TempZoom=False):

    for ind,t in enumerate(timeseries):
        print(t)
        data = dataset.sel(time=[t])
        temp = dataLVL2['ta'].sel(time=[t])
        #- select Tvals for the right axis:
        if TempZoom == True:
            tempArray = [0,-5,-10,-15,-20,25,-30,-35]
        else:
            tempArray = [0,-5,-10,-15,-20,-30,-40,-50,-60]
        TaxisVals,TaxisRange = get_closest_T(temp,tempArray)
        #- plot the data
        fig, axes = plt.subplots(ncols=3,figsize=(30,10))
        radData = {'XSpecH':{'data':10*np.log10(data['XSpecH']), 'axis':axes[0], 'lim':(-40,10)},
                   'KaSpecH':{'data':10*np.log10(data['KaSpecH']), 'axis':axes[1], 'lim':(-40,10)},
                   'WSpecH':{'data':10*np.log10(data['WSpecH']), 'axis':axes[2], 'lim':(-40,10)}}
        for rad in radData.keys():
            plot = radData[rad]['data'].plot(ax=radData[rad]['axis'],
                                        vmax=radData[rad]['lim'][1],
                                        vmin=radData[rad]['lim'][0],
                                        cmap='nipy_spectral',
                                        cbar_kwargs={'label':r'[dB]','pad':0.15})            
            
            secax = set_secax(radData[rad]['axis'])
            radData[rad]['axis'].set_yticks(TaxisRange)
            radData[rad]['axis'].set_yticklabels(tempArray)
            time = pd.to_datetime(str(t)).strftime('%Y%m%d %H:%M:%S')
            radData[rad]['axis'].set_title(rad +' '+ time, fontsize=14)
            plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
            radData[rad]['axis'].grid()
            radData[rad]['axis'].set_xlabel(r'Doppler velocity [ms$^{-1}]$',fontsize=14)
            radData[rad]['axis'].set_ylabel('Temp [K]',fontsize=14)
            radData[rad]['axis'].tick_params(axis='both', which='major', labelsize=12)
            if TempZoom == True:
                radData[rad]['axis'].set_ylim((0,TaxisRange[-1]))
                plotID = plotId+'_tempzoom'
                radData[rad]['axis'].set_xlim((-3.2,1.2))
                folder='zoom'
            else:
                plotID = plotId
                radData[rad]['axis'].set_xlim((-7,3))
                folder=''
        timePlotName = pd.to_datetime(str(t)).strftime('%Y%m%d_%H%M%S')
        plotFileName = ('_').join([timePlotName,plotID])
        year = pd.to_datetime(str(t)).strftime('%Y')
        month =  pd.to_datetime(str(t)).strftime('%m')
        day =  pd.to_datetime(str(t)).strftime('%d')
        filePathName = ('/').join([plotOutPath,year,month,day,folder,plotFileName])
        plt.savefig(filePathName+'.png',dpi=200, bbox_inches='tight')
        plt.close()
def get_closest_T(Tdata,T2find):
    Trange = np.zeros(len(T2find))
    Tvals  = np.zeros(len(T2find)) 
    for ind,temp in enumerate(T2find):
        Tind = (np.array((np.abs(Tdata.values-temp)).argmin()))
        Tvals[ind] = str(Tdata[0,Tind].values)
        Trange[ind] =  Tdata[0,Tind].range.values
    return Tvals,Trange
def set_secax(axis):
    def xrange(x):
       return x
    secax = axis.secondary_yaxis('right',functions=(xrange, xrange))
    secax.set_ylabel('height [m]', fontsize=14)
    secax.tick_params(axis='y', which='major', labelsize=12)
    return secax


def plotPol(data,plotOutPath, strDate, plotID, colmap='gist_ncar'):
    fig, axes = plt.subplots(nrows=4, figsize=(18,12))
    radData = {'ZDR':{'data':data, 'axis':axes[0], 'lim':(-0.5,3),'cmap':colmap},
               'DBZ':{'data':data, 'axis':axes[1], 'lim':(-30,15),'cmap':colmap},
               'KDP':{'data':data, 'axis':axes[2], 'lim':(-1,4), 'cmap':colmap},
               'RHV':{'data':data, 'axis':axes[3], 'lim':(0.97,1.001), 'cmap':colmap}}
    for rad in radData.keys():
        plot = radData[rad]['data'][rad].T.plot(ax=radData[rad]['axis'],
                                           vmax=radData[rad]['lim'][1],
                                           vmin=radData[rad]['lim'][0],
                                           cmap=radData[rad]['cmap'])

        radData[rad]['axis'].set_title(rad +' '+ strDate)
        plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
        radData[rad]['axis'].grid()
        radData[rad]['axis'].set_xlabel('')
        plt.tight_layout()
    plotFileName = ('_').join([strDate,plotID])
    filePathName = ('/').join([plotOutPath,plotFileName])
    plt.savefig(filePathName+'.png',dpi=200, bbox_inches='tight')

def testMo(data,plotOutPath, strDate, plotID, colmap='gist_ncar'):
    fig, axes = plt.subplots(nrows=3, figsize=(18,12))
    radData = {'ZDR':{'data':data, 'axis':axes[0], 'lim':(-0.5,3),'cmap':colmap},
               'ZH':{'data':data, 'axis':axes[1], 'lim':(-30,15),'cmap':'jet'},
               'ZV':{'data':data, 'axis':axes[2], 'lim':(-30,15), 'cmap':'jet'}}
    for rad in radData.keys():
        plot = radData[rad]['data'][rad].T.plot(ax=radData[rad]['axis'],
                                           vmax=radData[rad]['lim'][1],
                                           vmin=radData[rad]['lim'][0],
                                           cmap=radData[rad]['cmap'])

        radData[rad]['axis'].set_title(rad +' '+ strDate)
        plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
        radData[rad]['axis'].grid()
        radData[rad]['axis'].set_xlabel('')
        plt.tight_layout()
    plotFileName = ('_').join([strDate,plotID])
    filePathName = ('/').join([plotOutPath,plotFileName])
    plt.savefig(filePathName+'.png',dpi=200, bbox_inches='tight')
    plt.show()
def plotPolSpectra(dataset,dataLVL2,timeseries, plotOutPath, plotID):

    for ind,t in enumerate(timeseries):
        print(t)
        data = dataset.sel(time=[t])
        data2 = dataLVL2.sel(time=[t]) 
        #- plot the data
#        fig,axes = plt.subplots(nrows=1,ncols=1,figsize=(18,12))
        plt.pcolormesh(data['Vel2ZeroH'].fillna(0).values[0],
                              data2['ta'].values[0],
                              data['sZDR'].values[0],
                              vmin=-0.5,vmax=3,cmap=getNewNipySpectral())
        plt.colorbar(label='[dB]')
        time = pd.to_datetime(str(t)).strftime('%Y%m%d %H:%M:%S')
        plt.title('sZDR '+ time, fontsize=14)
        plt.grid()
        plt.ylabel(r'T [$^{\circ}$C]')
        plt.xlabel(r'Doppler velocity [ms$^{-1}]$',fontsize=14)
        #axes.tick_params(axis='both', which='major', labelsize=12)
        #plt.ylim([0,5000])
        plt.ylim([-5,-50])
        plt.xlim([-10,1])
        timePlotName = pd.to_datetime(str(t)).strftime('%Y%m%d_%H%M%S')
        plotFileName = ('_').join([timePlotName,plotID])
        year = pd.to_datetime(str(t)).strftime('%Y')
        month =  pd.to_datetime(str(t)).strftime('%m')
        day =  pd.to_datetime(str(t)).strftime('%d')
        filePathName = ('/').join([plotOutPath,year,month,day,plotFileName])
        plt.savefig(filePathName+'.png',dpi=200, bbox_inches='tight')
        #plt.show()
        plt.close()
 
def getNewNipySpectral():

    from matplotlib import cm
    from matplotlib.colors import ListedColormap

    numEnt = 15

    viridis = cm.get_cmap('nipy_spectral', 256)
    newcolors = viridis(np.linspace(0, 1, 256))

    colorSpace = np.linspace(198, 144, numEnt)/256
    colorTest=np.zeros((numEnt,4))
    colorTest[:,3] = 1
    colorTest[:,0]=colorSpace

    newcolors[- numEnt:, :] = colorTest
    newcmp = ListedColormap(newcolors)

    return newcmp

