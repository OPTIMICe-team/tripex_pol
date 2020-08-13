#####################################################
# This script is meant to postprocess the polarimetric spectra that were measured during the tripex-pol Campaign. For now it calculates the maximum sZDR and the 90% edges of the highest ZDR
#####################################################

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import glob as glob
from sys import argv
from netCDF4 import Dataset
from matplotlib import pyplot as plt
#own routines:
import processData as pro

#debugging==True

dateStart = pd.to_datetime('20190122'); dateEnd = pd.to_datetime('20190131')
dateList = pd.date_range(dateStart, dateEnd,freq='D')
dataOutPath = '/work/lvonterz/tripex_pol/output/'

for dayIndex,date2proc in enumerate(dateList):
    date2end = dateList[dayIndex+1]-pd.offsets.Second(4)
    timeFreq = '1H'; timeHour = pd.date_range(date2proc,date2end,freq=timeFreq)
    dataPath = dataOutPath + '/{year}/{month}/{day}'.format(year=date2proc.strftime('%Y'),
                                                            month=date2proc.strftime('%m'),
                                                            day=date2proc.strftime('%d'))
    for hourInd, hour in enumerate(timeHour[15:-1]):
        fileId = '{dateStr}_{HH}_tripex_pol_poldata_L0_spec_regridded.nc'.format(dateStr=hour.strftime('%Y%m%d'),HH=hour.strftime('%H'))
        file2proc = sorted(glob.glob(os.path.join(dataPath, fileId)))
        print((os.path.join(dataPath, fileId)))
        print(file2proc) 
        if file2proc:
            print('now processing: ',hour)
 #           if debugging==True:
 #               print(file2proc)
            data = xr.open_dataset(file2proc[0])
            #-- first mask according to snr:
            data = data.where(data.sSNR_H > 10.0)
            #data = data.where(data.sSNR_V > 10.0)
            #-- now getting the maximum sZDR along the velocity dimension:
            maxZDR = data['sZDR'].max(dim='Vel',keep_attrs=True)
            #-- now also get the index of the max: argmax itself has a problem with encountering all-NaN slices, so I need to get time-height of all nan velocitys:
            allna = data.sZDR.isnull().all(dim='Vel')
            idxMax = data.sZDR.where(~allna, 0).argmax(dim='Vel') # temporarely fill the allNans with 0 and get the idxmaxs
            idxMax = idxMax.where(~allna, np.nan) # perhaps instead of 0 you can fill again with NaNs the time-height where you had all nans
            #-- now we want to have the velocities where 10%(maxZDR) is reached (right and left of the maxZDR), to see how wide the enhanced ZDR are.
            maxZDR_m10 = maxZDR - 0.25*maxZDR
            # since we want to have it left and right of maximum, I am going to mask the right side of maximum when looking for the index on the left side and vice versa:
            dataMax = data.Vel[data.sZDR.where(~allna, 0).argmax(dim='Vel')]
            print(dataMax)
            maskedRight = data.sZDR.where(data.Vel>dataMax, np.nan)
            allna_right = maskedRight.isnull().all(dim='Vel')
            #maskedRight.isel(time=[1]).plot(cmap='nipy_spectral',vmin=-0.5,vmax=3)
            maskedLeft = data.sZDR.where(data.Vel<dataMax, np.nan)
            allna_left = maskedLeft.isnull().all(dim='Vel')
            # now calculate the left and right borders: 
            maxZDR_m10_diffs_right = np.abs(maskedRight-maxZDR_m10)
            idxMax_m10_right = maxZDR_m10_diffs_right.where(~allna_right,9999).argmin(dim='Vel')
            idxMax_m10_right = idxMax_m10_right.where(~allna_right, np.nan)
            
            maxZDR_m10_diffs_left = np.abs(maskedLeft-maxZDR_m10)
            idxMax_m10_left = maxZDR_m10_diffs_left.where(~allna_left,9999).argmin(dim='Vel')
            idxMax_m10_left = idxMax_m10_left.where(~allna_left, np.nan)

            dataMax_m10_left = data.Vel[maxZDR_m10_diffs_left.where(~allna_left,9999).argmin(dim='Vel')]
            dataMax_m10_right = data.Vel[maxZDR_m10_diffs_right.where(~allna_right,9999).argmin(dim='Vel')]
            dataMax_m10_left = dataMax_m10_left.where(~allna_left, np.nan)
            dataMax_m10_right = dataMax_m10_right.where(~allna_right, np.nan)
            print('idxMax_left',idxMax_m10_left)
            print('idxMax',idxMax) 
            print('idxMax_right',idxMax_m10_right)
            ##################################################################
            #-- now lets try a different method: maximum of second derivative
            ##################################################################
            dv = 0.01
            newVel = np.arange(np.min(data.Vel), np.max(data.Vel), dv)
           
            for t in data.time.values[1:-1]:
                VSpecInt = data.VSpec.sel(time=[t]).interp(Vel=newVel)
                VSpecInt = VSpecInt#*0.01
                HSpecInt = data.HSpec.sel(time=[t]).interp(Vel=newVel)
                HSpecInt = HSpecInt#*0.01
                sZDRInt = 10*np.log10(HSpecInt) - 10*np.log10(VSpecInt)
                sZDRsmooth = sZDRInt.rolling(Vel=10,min_periods=1,center=True).mean()
                firstDiff = sZDRsmooth.differentiate('Vel')
                firstDiffsmooth = firstDiff.rolling(Vel=10,min_periods=1,center=True).mean()
                secondDiff = firstDiffsmooth.differentiate('Vel')
               
                #- now lets calculate maximum of 1st derivative (this should be the "wendepunkt" of the maximum of sZDR)
                maxFirstDiff = firstDiffsmooth.max(dim='Vel')
                allna_Diff = firstDiffsmooth.isnull().all(dim='Vel')
                idxMaxFirstDiff = firstDiffsmooth.where(~allna_Diff, 0).argmax(dim='Vel') # temporarely fill the allNans with 0 and get the idxmaxs
                idxMaxFirstDiff = idxMaxFirstDiff.where(~allna_Diff, np.nan) # perhaps instead of 0 you can fill again with NaNs the time-height where you had all nans
                idxMaxFirstDiffVel = firstDiffsmooth.Vel[firstDiffsmooth.where(~allna_Diff, 0).argmax(dim='Vel')]
                idxMaxFirstDiffVel = idxMaxFirstDiffVel.where(~allna_Diff, np.nan)
                print(maxFirstDiff)
                #-- now maximum of second derivatice (which should be the onset of the maximum of sZDR?)
                maxsecondDiff = secondDiff.max(dim='Vel')
                allna_Diff_second = secondDiff.isnull().all(dim='Vel')
                idxMaxSecondDiff = secondDiff.where(~allna_Diff_second, 0).argmax(dim='Vel') # temporarely fill the allNans with 0 and get the idxmaxs
                idxMaxSecondDiff = idxMaxSecondDiff.where(~allna_Diff_second, np.nan) # perhaps instead of 0 you can fill again with NaNs the time-height where you had all nans
                idxMaxSecondDiffVel = secondDiff.Vel[secondDiff.where(~allna_Diff_second, 0).argmax(dim='Vel')]
                idxMaxSecondDiffVel = idxMaxSecondDiffVel.where(~allna_Diff_second, np.nan)
                
                data.sZDR.isel(time=[1],height=[50]).plot()
                plt.axvline(x=dataMax.isel(time=[1],height=[50]).values[0],color='k',label='max')
                plt.axvline(x=dataMax_m10_left.isel(time=[1],height=[50]).values[0],color='grey',label='-25%')
                plt.axvline(x=dataMax_m10_right.isel(time=[1],height=[50]).values[0],color='grey',label='+25%')
                plt.axvline(x=idxMaxFirstDiffVel.isel(height=[50]).values[0],color='g',label='max 1st deriv')
                plt.axvline(x=idxMaxSecondDiffVel.isel(height=[50]).values[0],color='r',label='max 2nd deriv')
                plt.legend()
                plt.savefig('maxInd_test_height50.png')
                plt.show()
                quit()



                #firstDiff = sZDRsmooth.diff(dim='Vel')/dv
                #secondDiff = firstDiff.diff(dim='Vel')/dv        
                #data.sZDR.isel(time=[1],height=[50]).plot(label='sZDR')
                #sZDRInt.isel(height=[50]).plot(label='int')
                #firstDiff.isel(height=[50]).plot(label='1st_diff')
                #secondDiff.isel(height=[50]).plot(label='2nd_diff')
                plt.legend()
                #plt.savefig('derivative_smoothing.png')
                plt.show()
                quit()

            data.sZDR.isel(time=[1],height=[50]).plot()
            plt.axvline(x=dataMax.isel(time=[1],height=[50]).values[0],color='k',label='max')
            plt.axvline(x=dataMax_m10_left.isel(time=[1],height=[50]).values[0],color='grey',label='-25%')
            plt.axvline(x=dataMax_m10_right.isel(time=[1],height=[50]).values[0],color='grey',label='+25%')
            plt.legend()
            plt.savefig('maxInd_test_height50.png')
            plt.show()
            quit()
#            data.sZDR.isel(time=[1]).plot(vmin=-0.5,vmax=3,cmap='nipy_spectral')
#            plt.plot(dataMax.isel(time=[1]).values[0],data.height.values,'k',linewidth=1.5)
#            plt.plot(dataMax_m10_left.isel(time=[1]).values[0],data.height.values,'k--',linewidth=1.5)
#            plt.plot(dataMax_m10_right.isel(time=[1]).values[0],data.height.values,'k--',linewidth=1.5)
#            plt.ylim([0,4000])
#            plt.xlim([-6,6])
#            plt.savefig('maxInd_test.png')
#            plt.show()
#            quit()
            
'''
            for t in data.time.values[1:-1]:
                datasel = data.sel(time=[t])
                maxZDR = datasel['sZDR'].max(dim='Vel',keep_attrs=True)
                maxZDRind = np.ones_like(maxZDR.values[0])*np.nan
                print(maxZDRind.shape)
                maxZDR_p10 = maxZDR + maxZDR/10
                maxZDR_m10 = maxZDR - maxZDR/10
                print(maxZDR_m10.values)
                maxZDR_m10 = np.abs(maxZDR.values-maxZDR_m10.values).min()
                print(maxZDR_m10)
                quit()
                for ind,r in enumerate(datasel.height.values[1:-1]):
                    dataselra = datasel.sel(height=[r]) 
                    maxInd = dataselra['sZDR'].where(dataselra.sZDR==maxZDR,drop=True).squeeze()
                    maxInd_p10 = dataselra['sZDR'].where(dataselra.sZDR==maxZDR_p10,drop=True).squeeze()
                    maxInd_m10 = dataselra['sZDR'].where(np.abs(dataselra.sZDR-maxZDR_m10).min(),drop=True).squeeze()
                    print(maxInd.Vel.values)
                    #maxInd = dataselra['sZDR'].where(dataselra.sZDR==dataselra.sZDR.max(dim='Vel'),drop=True).squeeze()
                    print(maxInd_p10.Vel.values)
                    print(maxInd_m10.Vel.values)
                    quit()
                    if maxInd:                  
                        maxZDRind[ind] = maxInd.Vel.values
                    else:
                        maxZDRind[ind] = np.nan
                #dataselDF = datasel['sZDR'].to_dataframe()
                #dataselDF=dataselDF.groupby(level=1).idxmax()
                #print(maxZDR)
                #dataselDFval=dataselDF.values[1]
                #print(dataselDFval.type)
                #maxZDRind = datasel['sZDR'].idxmax(dim='Vel',skipna=True)
                datasel['sZDR'].plot(cmap='nipy_spectral',vmin=-0.5,vmax=3)
                plt.plot(maxZDRind,datasel.height.values,'k',linewidth=1)
                plt.ylim([0,4000])
                plt.savefig('maxInd_test.png')
                plt.show()
                quit()
'''
