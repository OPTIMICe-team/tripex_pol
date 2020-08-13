import glob
import numpy as np
import pandas as pd
from sys import argv
from netCDF4 import Dataset

import writeData
import tripexLib as trLib
import readRadarInfo as rdInfo


#input File Path
path = argv[1]
#output File Path
outputPath = argv[2]
#prefix
prefixL1 = argv[3]
#-----------------------------

#--Time Definitions-----------
yearStr = argv[4]
year = int(yearStr)
monthStr = argv[5]
month = int(monthStr)
dayStr = argv[6]
day = int(dayStr)
beguinTimeRefStr = argv[7] 
beguinTimeRef = int(beguinTimeRefStr)
timeFreq = argv[8]
timeTolerance = argv[9]

dateName = yearStr + monthStr + dayStr
endTimeRef = beguinTimeRef + 1
start = pd.datetime(year, month, day,beguinTimeRef, 0, 0)

if endTimeRef == 24:
   end = pd.datetime(year, month, day, 23, 59, 59)
else:
   end = pd.datetime(year, month, day, endTimeRef, 0, 0)

timeRef = pd.date_range(start, end, freq=timeFreq)
timeRefUnix = np.array(timeRef,float)
usedIndexTime = np.ones((len(timeRef)))*np.nan
#-----------------------------

#--Range Definitions----------
beguinRangeRef = int(argv[10])
endRangeRef = int(argv[11])
rangeFreq = int(argv[12])
rangeTolerance = int(argv[13])
rangeRef = np.arange(beguinRangeRef, endRangeRef, rangeFreq)
usedIndexRange = np.ones((len(rangeRef)))*np.nan
#-----------------------------

#--Radar variables------------
radar = argv[14]
rangeGateOffSet = float(argv[15]) 
variableName = argv[16] 
#zeOffset = float(argv[17])
#-----------------------------

#output File Definitions
outPutFile = ('_').join([prefixL1, dateName,
                       beguinTimeRefStr+'.nc'])
outPutFilePath = ('/').join([outputPath, outPutFile])
print outPutFilePath

#output variable name
if variableName == 'Zg':    
   varFinalName = 'Ze'
 
elif variableName == 'vm' or variableName == 'vd':
   varFinalName = 'v'

elif variableName == 'VELg':
   varFinalName = 'v'

elif variableName == 'RMS':
   varFinalName = 'SW'

elif variableName == 'sigma':
   varFinalName = 'SW'

elif variableName == 'W':
   varFinalName = 'SW'

else:
   varFinalName = variableName

varNameOutput = ('_').join([varFinalName, radar])

#--------------------------------------------
#Files to process
fileList = rdInfo.getFileList(radar, year, month,
                             day, beguinTimeRef)
#fileList = trLib.getFileList(path, dateName,
#			    beguinTimeRef, radar)
#
#if radar == 'Ka':
#   fileList = trLib.checkFileListKa(fileList)
#--------------------------------------------


varData = pd.DataFrame()
for radarFile in fileList:

   #it opens the file 
   print radarFile
   try:
      rootgrp = Dataset(radarFile, 'r')
      
   except:
      print 'No Data'
      break
   epoch = trLib.getEpochTime(rootgrp, radar)
   timesW = rootgrp.variables['time'][:]
  
   #timesW = timesWAtt[:]
   humamTimeW = epoch + pd.to_timedelta(timesW, unit='S')
   print humamTimeW 
   quit()
   #it gets the range and corrects for the reference height
   ranges =rootgrp.variables['range'][:] + rangeGateOffSet  
    
   #it gets desireble variable 
   var = rootgrp.variables[variableName][:]

   #Ka version ------------------------------
   if varFinalName == 'Ze' and radar == 'W':
      var = 10*np.log10(var)
      print varFinalName, radar
   else:
      pass

   if varFinalName == 'Ze' and radar == 'Ka':
      var = 10*np.log10(var)
   else:
      pass

   if varFinalName == 'LDRg' and radar == 'Ka':
      var = 10*np.log10(var)
   else:
      pass
   #-----------------------------------------

   if radar == 'W':
      var[var==-999.] = np.nan 

   if radar == 'X' and varFinalName == 'Ze':
      var[var==-32] = np.nan
 #  if radar == 'X' and varFinalName == 'ZDR':
 #     var[var==-32] = np.nan
 #  if radar == 'X' and varFinalName == 'KDP':
 #     var[var==-32] = np.nan
 #  if radar == 'X' and varFinalName == 'PhiDP':
 #     var[var==-32] = np.nan
 #  if radar == 'X' and varFinalName == 'RhoHV':
 #     var[var==-32] = np.nan


   var = np.ma.masked_invalid(var)
   np.ma.set_fill_value(var, np.nan) 
  
   varDataTemp = pd.DataFrame(index=humamTimeW, columns=ranges, data=var)
   varDataTemp['times'] = humamTimeW 
   varDataTemp = varDataTemp.sort_values(by=['times'], ascending=[True])
   varDataTemp = varDataTemp.drop_duplicates(subset=['times'])

   varData=varData.append(varDataTemp) 
   rootgrp.close()

try:
   #resample in time
   timeIndexList = trLib.getIndexList(varData, timeRef, timeTolerance)
   emptyDataFrame = pd.DataFrame(index=timeRef, columns=varData.columns)
   resampledTime = trLib.getResampledDataPd(emptyDataFrame, varData, timeIndexList)
   timeDeviation = trLib.getDeviationPd(timeRef, resampledTime, timeTolerance)
   del resampledTime['times']
   resampledTime = resampledTime.T

except:
   print 'No Times'

try:
   #resample in range
   resampledTime['ranges']=ranges
   rangeIndexList = trLib.getIndexList(resampledTime, rangeRef, rangeTolerance)
   emptyDataFrame = pd.DataFrame(index=rangeRef, columns=resampledTime.columns)
   resampledTimeRange = trLib.getResampledDataPd(emptyDataFrame, resampledTime,
                                             rangeIndexList) 
   rangeDeviation = trLib.getDeviationPd(rangeRef, resampledTimeRange)
   del resampledTimeRange['ranges'] 

except:
   print 'No Ranges'

timeRefUnix = timeRefUnix/10.**9

try:  
   rootgrpOut = writeData.createNetCdf(outPutFilePath, prefixL1)
   time_ref = writeData.createTimeDimension(rootgrpOut, timeRefUnix, prefixL1)
   range_ref = writeData.createRangeDimension(rootgrpOut, rangeRef, prefixL1)

   timeDeviation = np.array(timeDeviation.astype(np.float32))
   time_dev = writeData.createDeviation(rootgrpOut, timeDeviation,
                                    'delta_time', radar, prefixL1)

   rangeDeviation = np.array(rangeDeviation.astype(np.float32))
   range_dev = writeData.createDeviation(rootgrpOut, rangeDeviation,
                                          'delta_altitude', radar, prefixL1)

   if varNameOutput == 'Ze_'+radar:
      resampledTimeRange = resampledTimeRange + zeOffset

   resampledTimeRange = np.array(resampledTimeRange.T.astype(np.float32))
   var_resampled = writeData.createVariable(rootgrpOut, resampledTimeRange,
                                        varFinalName, varNameOutput, radar, 
                                        prefixL1, np.float32)
except:
#except IOError:
   print 'No Data To Write :('

rootgrpOut.close()

print 'Done %s'%(varNameOutput)
