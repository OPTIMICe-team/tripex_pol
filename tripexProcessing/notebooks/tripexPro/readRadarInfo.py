import pandas as pd
import numpy as np
import datetime

def getFileList(radar, year, month, day, beguinTimeRef):
   outputPath = '/home/lvonterz/tripexPro/outputtest' #'/home/jdias/Projects/radarDataResampled/radarInfo'

   #input
   #radar = 'X'
   #year = 2015
   #month = 11
   #day = 24
   #beguinTimeRef = 17

   endTimeRef = beguinTimeRef + 1
   start = pd.datetime(year, month, day,beguinTimeRef, 0, 0)

   if endTimeRef == 24: 
      end = pd.datetime(year, month, day, 23, 59, 59) 
   else:
      end = pd.datetime(year, month, day, endTimeRef, 0, 0)

   if radar == 'Ka':
      outputInfoName = 'kaBandInfo.csv'
      outputFile = ('/').join([outputPath, outputInfoName])
      readData=pd.read_csv(outputFile, index_col='Unnamed: 0')
      date=readData.time_long_name.str.split(' ').str[-3:]#.str.split('.')
    
   if radar == 'W':
      outputInfoName = 'wBandInfo.csv'
      outputFile = ('/').join([outputPath, outputInfoName])
      readData=pd.read_csv(outputFile, index_col='Unnamed: 0')    
      date=readData.time_long_name.str.split(' ').str[-2:]
    
   if radar == 'X':
      outputInfoName = 'xBandInfo.csv'
      outputFile = ('/').join([outputPath, outputInfoName])
      readData=pd.read_csv(outputFile, index_col='Unnamed: 0')
      date=readData.time_long_name.str.split(' ').str[-2:]

   epoch = pd.to_datetime(date.str.get(0))
   readData = readData.sort_values(by=['filePath'], ascending=[True])
   humanDateMin=pd.to_timedelta(readData.min_time, unit='s')+epoch
   humanDateMax=pd.to_timedelta(readData.max_time, unit='s')+epoch
#  humanDateMin=pd.to_datetime(humanDateMin)
#  humanDateMax=pd.to_datetime(humanDateMax)
   readData['humanDateMin']=humanDateMin
   readData['humanDateMax']=humanDateMax
   readData['year']=humanDateMin.dt.year
   readData['month']=humanDateMin.dt.month
   readData['day']=humanDateMin.dt.day
   readData['hourMin'] = humanDateMin.dt.hour
   readData['hourMax'] = humanDateMax.dt.hour

   if radar == 'X':
      #fazer a busca com tempo em segundos
      #X band
      tempFrame = readData[readData.elv==90.]
      tempFrame = tempFrame[tempFrame.year==year]
      tempFrame = tempFrame[tempFrame.month==month]
      tempFrame = tempFrame[tempFrame.day==day]

   elif radar == 'W':
     #fazer a busca com tempo em segundos
     #W band
     tempFrame = readData[readData.elv==90.]
     tempFrame = tempFrame[tempFrame.humanDateMin>=start]
     tempFrame = tempFrame[tempFrame.humanDateMax<end]
     tempFrame = tempFrame[tempFrame.day==day]
     tempFrame = tempFrame[tempFrame.hourMin==beguinTimeRef]


   else:
     tempFrame = readData[readData.elv==90.]
     tempFrame = tempFrame[tempFrame.year==year]
     tempFrame = tempFrame[tempFrame.month==month]
     tempFrame = tempFrame[tempFrame.day==day]
     tempFrame = tempFrame[tempFrame.hourMin==beguinTimeRef]

     if len(tempFrame)==0:
        tempFrame = readData[readData.elv==90.]
        tempFrame = tempFrame[tempFrame.year==year]
        tempFrame = tempFrame[tempFrame.month==month]
        tempFrame = tempFrame[tempFrame.day==day]
        tempFrame = tempFrame[tempFrame.hourMin<=beguinTimeRef]
        tempFrame = tempFrame[tempFrame.hourMax>=beguinTimeRef]


   finalFileList = tempFrame.filePath.values
   return list(finalFileList)

