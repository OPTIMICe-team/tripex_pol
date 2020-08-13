import numpy as np
import pandas as pd
import glob
from netCDF4 import Dataset

def getEmptyMatrix(rows, cols):

   emptyMatrix = np.ones((rows, cols))
   emptyMatrix[:] = np.nan
   emptyMatrix = np.ma.masked_invalid(emptyMatrix)
    
   return emptyMatrix

def getEpochTime(rootgrp, radar):

   if radar == 'W':
      long_name = getattr(rootgrp.variables['time'], 'long_name').split(' ')
      
      date = long_name[-2]
      time = long_name[-1]
            
   elif radar == 'X':
      long_name = getattr(rootgrp.variables['time'], 'long_name').split(' ')  
      date = long_name[-2]
      time = long_name[-1]

   else:
      long_name = getattr(rootgrp.variables['time'], 'long_name').split(' ')
      date = long_name[-3]
      time = long_name[-2]
         
   epoch = pd.to_datetime(' '.join([date, time]))
    
   return epoch

   
def getFileList(dataPath, date, beguinTime, radar):
   
   
   beguinTime = str(beguinTime)
   if radar == 'W':
      prefix = 'joyrad94_joyce_compact_'
      suffix = ('').join([beguinTime,'.nc'])
      fileNameSearch = ('*').join([date, suffix])
      fileNameSearch = ('').join([prefix, fileNameSearch])
      fileNameSearch = ('/').join([dataPath, fileNameSearch])
    
   elif radar == 'X':
      suffix = 'kixpol.nc'
      fileNameSearch = ('*').join([date, suffix])
      fileNameSearch = ('/').join([dataPath, fileNameSearch])
        
   else:
      suffix = 'mmclx.00'
      fileNameSearch = ('*').join([beguinTime, suffix])
      fileNameSearch = ('_').join([date,fileNameSearch])
      fileNameSearch = ('/').join([dataPath, fileNameSearch])

   fileList = glob.glob(fileNameSearch)
   #fileList = fileList.sort()
   #print fileNameSearch
   return fileList

 

def getIndexList(dataTable, reference, tolerance):
    
   indexList = []
   for value in reference:
    
      index = getNearestIndex(dataTable, value, tolerance)
      indexList.append(index)

   return indexList


def getNearestIndex(timeRef, timeStamp, tolerance):

   try:
      index = timeRef.index.get_loc(timeStamp, method='nearest', 
                                   tolerance=tolerance)
    
   except:
      index = np.nan
    
   return index

def getDeviation(rangeRef, ranges, usedIndexRange):
   
   deviationList = []
   for i, usedIndex in enumerate(usedIndexRange):

      try:
         if type(rangeRef) == pd.tseries.index.DatetimeIndex:
            deviation = rangeRef[i].second - ranges[int(usedIndex)].second
         else:
            deviation = rangeRef[i] - ranges[int(usedIndex)]

         deviationList.append(deviation)
         
      except:
         deviationList.append(np.nan)
    
   return np.array(deviationList) 

def getResampledData(emptyMatrixData, variableData, indexList, usedIndex):
    
   usedIndex = []
   for i,index in enumerate(indexList):

      try:
         emptyMatrixData[i] = variableData[int(index)]
         usedIndex.append(int(index))
       
      except:
         usedIndex.append(np.nan)
    
   return emptyMatrixData, usedIndex

def getResampledDataPd(emptyDataFrame, dataDataFrame, indexList):

   for i, index in enumerate(indexList):

      try:
         emptyDataFrame.iloc[i]=dataDataFrame.iloc[index]
        
      except:
         pass
        
   return emptyDataFrame

def getDeviationPd(reference, resampledData, tolerance=None):

   if type(reference) == pd.DatetimeIndex:

      deviation = reference.second-pd.to_datetime(list(resampledData.times)).second
      resampledData['delta_time']=deviation
      tolerance = int(tolerance[:-1])
      correction = 60.+resampledData.delta_time[resampledData.delta_time<-tolerance]
      resampledData.loc[resampledData.delta_time<-tolerance,'delta_time']=correction
      deviation = resampledData.delta_time
      del resampledData['delta_time']     

   else:
      deviation = reference - resampledData.ranges

   return deviation


def checkFileListKa(fileList):

   listAux = []
   for nameFile in fileList:
      rootgrp = Dataset(nameFile, 'r')
      elv = rootgrp.variables['elv'][:]
      lenRange = len(rootgrp.variables['range'])
      if len(np.argwhere(elv !=90))==0:
         listAux.append(nameFile)
      rootgrp.close()

   return listAux

