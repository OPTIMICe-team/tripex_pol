##########################################################################################
# calculates and plots the moments (written mainly by Jose, I have my own plotting routine)
##########################################################################################

import os
import time 
import xarray as xr
import numpy as np 
import netCDF4 as nc
import glob as gb

import matplotlib.pyplot as plt

from sys import argv

startTime = time.time()

def convertToXrDtArr(variable, times, ranges):

    varXr = xr.DataArray(variable, 
                         dims=('time', 'range'),
                         coords={'time':times, 'range':ranges})

    return varXr


def globalAttributes(data):
    
    data.attrs['Experiment']= 'TRIPEX-POL, Forschungszentrum Juelich'
    data.attrs['Instrument']= 'GRARAD-94 94.4 GHz cloud radar (W band)'
    data.attrs['Data']= 'Produced by Jose Dias, jdiasnet@uni-koeln.de'
    data.attrs['Institution']= 'Data processed within the Emmy-Noether Group OPTIMIce, Institute for Geophysics and Meteorology, University of Cologne, Germany'
    data.attrs['Latitude']= '50.908547 N'
    data.attrs['Longitude']= '6.413536 E'
    data.attrs['Altitude']= 'Altitude of the JOYCE (www.joyce.cloud) platform: 111m asl'
    
    return data


def variableAtrributes(data, band):

    data.range.attrs['long_name']= 'Range from antenna to the center of each range gate'
    data.range.attrs['units']= 'm'
    
    data.time.attrs['long_name']='time'
    data.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
    if band == 'X':
        data.X_DBZ_V.attrs['standard_name']= 'equivalent_reflectivity_factor'
        data.X_DBZ_V.attrs['long_name']= 'equivalent reflectivity factor, '+band+' Band'
        data.X_DBZ_V.attrs['units']= 'dBZ'
        data.X_VEL_V.attrs['standard_name']= 'radial_velocity_of_scatterers_away_from_instrument'
        data.X_VEL_V.attrs['long_name']= 'mean Doppler velocity, '+band+' Band'
        data.X_VEL_V.attrs['units']= 'm s-1'
        data.X_WIDTH_V.attrs['standard_name']= 'doppler_spectrum_width'
        data.X_WIDTH_V.attrs['long_name']= 'doppler spectrum width, '+band+' Band'
        data.X_WIDTH_V.attrs['units']= 'm s-1'
        data.X_SK_V.attrs['long_name']= 'skewness, '+band+' Band'
        data.X_SK_V.attrs['units']= '-'

    if band =='Ka':
        data.Ka_DBZ_H.attrs['standard_name']= 'equivalent_reflectivity_factor'
        data.Ka_DBZ_H.attrs['long_name']= 'equivalent reflectivity factor, horizontal polarization, '+band+' Band'
        data.Ka_DBZ_H.attrs['units']= 'dBZ'
        data.Ka_VEL_H.attrs['standard_name']= 'radial_velocity_of_scatterers_away_from_instrument'
        data.Ka_VEL_H.attrs['long_name']= 'mean Doppler velocity, horizontal polarization, '+band+' Band'
        data.Ka_VEL_H.attrs['units']= 'm s-1'
    	
        data.Ka_DBZ_V.attrs['standard_name']= 'equivalent_reflectivity_factor'
        data.Ka_DBZ_V.attrs['long_name']= 'equivalent reflectivity factor, vertical polarization, '+band+' Band'
        data.Ka_DBZ_V.attrs['units']= 'dBZ'
        data.Ka_VEL_V.attrs['standard_name']= 'radial_velocity_of_scatterers_away_from_instrument'
        data.Ka_VEL_V.attrs['long_name']= 'mean Doppler velocity, vertical polarization, '+band+' Band'
        data.Ka_VEL_V.attrs['units']= 'm s-1'
        
        data.Ka_WIDTH_V.attrs['standard_name']= 'doppler_spectrum_width'
        data.Ka_WIDTH_V.attrs['long_name']= 'radar spectral width, vertical polarization'
        data.Ka_WIDTH_V.attrs['units']= 'm s-1'
        data.Ka_WIDTH_H.attrs['standard_name']= 'doppler_spectrum_width'
        data.Ka_WIDTH_H.attrs['long_name']= 'radar spectral width, horizontal polarization'
        data.Ka_WIDTH_H.attrs['units']= 'm s-1'
        
        data.Ka_SK_H.attrs['long_name']= 'skewness, horizontal polarization, X Band'
        data.Ka_SK_H.attrs['units']= '-'
        data.Ka_SK_V.attrs['long_name']= 'skewness, vertical polarization'
        data.Ka_SK_V.attrs['units']= '-'
    
    if band =='W':
        data.W_DBZ_V.attrs['standard_name']= 'equivalent_reflectivity_factor'
        data.W_DBZ_H.attrs['long_name']= 'equivalent reflectivity factor, horizontal polarization, '+band+' Band'
        data.W_DBZ_H.attrs['units']= 'dBZ'
        data.W_VEL_H.attrs['standard_name']= 'radial_velocity_of_scatterers_away_from_instrument'
        data.W_VEL_H.attrs['long_name']= 'mean Doppler velocity, horizontal polarization, '+band+' Band'
        data.W_VEL_H.attrs['units']= 'm s-1'
    	
        data.W_DBZ_H.attrs['standard_name']= 'equivalent_reflectivity_factor'
        data.W_DBZ_V.attrs['long_name']= 'equivalent reflectivity factor, vertical polarization, '+band+' Band'
        data.W_DBZ_V.attrs['units']= 'dBZ'
        data.W_VEL_V.attrs['standard_name']= 'radial_velocity_of_scatterers_away_from_instrument'
        data.W_VEL_V.attrs['long_name']= 'mean Doppler velocity, vertical polarization, '+band+' Band'
        data.W_VEL_V.attrs['units']= 'm s-1'
    
        data.W_WIDTH_V.attrs['standard_name']= 'doppler_spectrum_width'
        data.W_WIDTH_V.attrs['long_name']= 'radar spectral width, vertical polarization, '+band+' Band'
        data.W_WIDTH_V.attrs['units']= 'm s-1'
        data.W_WIDTH_H.attrs['standard_name']= 'doppler_spectrum_width'
        data.W_WIDTH_H.attrs['long_name']= 'radar spectral width, horizontal polarization, '+band+' Band'
        data.W_WIDTH_H.attrs['units']= 'm s-1'
        
        data.W_SK_H.attrs['standard_name']= 'doppler_spectrum_width'
        data.W_SK_H.attrs['long_name']= 'skewness, horizontal polarization, '+band+' Band'
        data.W_SK_H.attrs['units']= '-'
        data.W_SK_V.attrs['standard_name']= 'doppler_spectrum_width'
        data.W_SK_V.attrs['long_name']= 'skewness, vertical polarization, '+band+' Band'
        data.W_SK_V.attrs['units']= '-'
    

    return data


def getMomentsW(data):
# original one from Jose, I have written my own!
    #chirpNumber = np.int(data.attrs['chirp_number'])
    chirpNumber=4
    for chirp in range(chirpNumber):
    
        chVel = 'C{chirp}Vel'.format(chirp=chirp+1)
        chSpecV = 'C{chirp}specV'.format(chirp=chirp+1)
        chZeV = 'C{chirp}ZeV'.format(chirp=chirp+1)
        chRg = 'C{chirp}Range'.format(chirp=chirp+1)
        chSpecH = 'C{chirp}specH'.format(chirp=chirp+1)
        chZeH = 'C{chirp}ZeH'.format(chirp=chirp+1)
        
        #Ze 
        dv=1
        dv = data[chVel].values[1] - data[chVel].values[0]
        ZeH = data[chSpecH].sum(dim=chVel)*dv
        ZeV = data[chSpecV].sum(dim=chVel)*dv
        
        #MDV
        mdvH = ((data[chSpecH]*data[chVel]).sum(dim=chVel)*dv)/ZeH
        mdvV = ((data[chSpecV]*data[chVel]).sum(dim=chVel)*dv)/ZeV
        
        #SW
        sqrDiffV = (data[chVel] - mdvV)**2
        swV = (((data[chSpecV]*sqrDiffV).sum(dim=chVel)*dv)/ZeV)
        swV = swV**(1/2.)
        sqrDiffH = (data[chVel] - mdvH)**2
        swH = (((data[chSpecH]*sqrDiffH).sum(dim=chVel)*dv)/ZeH)
        swH = swH**(1/2.)

        #SK
        cubDiffV = (data[chVel] - mdvV)**3
        skV = (((data[chSpecV]*cubDiffV).sum(dim=chVel)*dv)/(ZeV*swV**3))
        cubDiffH = (data[chVel] - mdvH)**3
        skH = (((data[chSpecH]*cubDiffH).sum(dim=chVel)*dv)/(ZeH*swH**3))

        if chirp == 0:
            tempZeV = ZeV.values
            tempRange = data[chRg].values
            tempMDVV = mdvV.values
            tempSWV = swV.values
            tempSKV = skV.values
            tempZeH = ZeH.values
            tempMDVH = mdvH.values
            tempSWH = swH.values
            tempSKH = skH.values

        else:
            tempZeV = np.concatenate((tempZeV, ZeV.values), axis=1)
            tempRange = np.concatenate((tempRange, data[chRg].values))
            tempMDVV = np.concatenate((tempMDVV, mdvV.values), axis=1)
            tempSWV = np.concatenate((tempSWV, swV.values), axis=1)
            tempSKV = np.concatenate((tempSKV, skV.values), axis=1)
            tempZeH = np.concatenate((tempZeH, ZeH.values), axis=1)
            tempMDVH = np.concatenate((tempMDVH, mdvH.values), axis=1)
            tempSWH = np.concatenate((tempSWH, swH.values), axis=1)
            tempSKH = np.concatenate((tempSKH, skH.values), axis=1)


    tempZeV = 10*np.log10(tempZeV)
    tempZeH = 10*np.log10(tempZeH)
    tempMome = xr.Dataset({'WZeH':convertToXrDtArr(tempZeH, data.time, tempRange),
                           'WrvH':convertToXrDtArr(tempMDVH, data.time, tempRange),
                           'WswH':convertToXrDtArr(tempSWH, data.time, tempRange),
                           'WskH':convertToXrDtArr(tempSKH, data.time, tempRange),
			   'WZeV':convertToXrDtArr(tempZeV, data.time, tempRange),
                           'WrvV':convertToXrDtArr(tempMDVV, data.time, tempRange),
                           'WswV':convertToXrDtArr(tempSWV, data.time, tempRange),
                           'WskV':convertToXrDtArr(tempSKV, data.time, tempRange)})
    
    return tempMome

def getMoments(data,band):
#- this calculates the moments directly from spectrum. We need to specify the band because X only has co-channel and W-band needs to be divided by dv, since we regridded the Doppler vel. coordinates.
    if band == 'X':
        ZeV = data[band+'SpecV'].sum(dim='doppler'+band)       
    #MDV
        mdvV = ((data[band+'SpecV']*(data['doppler'+band])).sum(dim='doppler'+band))/ZeV
        
        sqrDiffV = (data['doppler'+band] - mdvV)**2
        swV = (((data[band+'SpecV']*sqrDiffV).sum(dim='doppler'+band))/ZeV)
        swV = swV**(1/2.)
        
        cubDiffV = (data['doppler'+band] - mdvV)**3
        skV = (((data[band+'SpecV']*cubDiffV).sum(dim='doppler'+band))/(ZeV*swV**3))
        ZeV = 10*np.log10(ZeV)
        tempMome = xr.Dataset({'X_DBZ_V':ZeV,
                           'X_VEL_V':mdvV,
                           'X_WIDTH_V':swV,
                           'X_SK_V':skV,
			   })
    else:
        if band == 'Ka':
            dv = 1
        elif band == 'W':
            dv = data['dvW'].values[0]
        ZeH = data[band+'SpecH'].sum(dim='doppler'+band)*dv
        ZeV = data[band+'SpecV'].sum(dim='doppler'+band)*dv
        mdvH = ((data[band+'SpecH']*(data['doppler'+band])).sum(dim='doppler'+band)*dv)/ZeH
        mdvV = ((data[band+'SpecV']*(data['doppler'+band])).sum(dim='doppler'+band)*dv)/ZeV
        
        sqrDiffV = (data['doppler'+band] - mdvV)**2
        swV = (((data[band+'SpecV']*sqrDiffV).sum(dim='doppler'+band)*dv)/ZeV)
        swV = swV**(1/2.)
        sqrDiffH = (data['doppler'+band] - mdvH)**2
        swH = (((data[band+'SpecH']*sqrDiffH).sum(dim='doppler'+band)*dv)/ZeH)
        swH = swH**(1/2.)
                
        cubDiffH = (data['doppler'+band] - mdvH)**3
        skH = (((data[band+'SpecH']*cubDiffH).sum(dim='doppler'+band)*dv)/(ZeH*swH**3))
        cubDiffV = (data['doppler'+band] - mdvV)**3
        skV = (((data[band+'SpecV']*cubDiffV).sum(dim='doppler'+band)*dv)/(ZeV*swV**3))
        ZeV = 10*np.log10(ZeV)
        ZeH = 10*np.log10(ZeH)
        tempMome = xr.Dataset({band+'_DBZ_H':ZeH,
                           band+'_VEL_H':mdvH,
                           band+'_WIDTH_H':swH,
                           band+'_SK_H':skH,
			   band+'_DBZ_V':ZeV,
                           band+'_VEL_V':mdvV,
                           band+'_WIDTH_V':swV,
                           band+'_SK_V':skV})
    
    return tempMome

def plotCleanMon(dataset, plotOutPath, strDate, band):
#-- plots moments (from Jose) I have my own in /net/broebroe/lvonterz/tripex_pol/plotting/...
    fig, axes = plt.subplots(nrows=3, figsize=(18,24))

    radData = {band+'ZeV':{'data':dataset, 'axis':axes[0], 'lim':(-35,25)},
               band+'rvV':{'data':dataset, 'axis':axes[1], 'lim':(-3,0)},
               band+'swV':{'data':dataset, 'axis':axes[2],'lim':(0,1)},
               #band+'skH':{'data':dataset, 'axis':axes[3],'lim':(-1,1)}
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

    plotFileName = ('_').join([strDate,'momTEST.png'])
    filePathName = ('/').join([plotOutPath,plotFileName])
    #plt.savefig(filePathName+'.png')
    #plt.close()
    #plt.show()

