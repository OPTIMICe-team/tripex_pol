import pandas as pd
import xarray as xr
import numpy as np
from netCDF4 import Dataset
import attenuationLib as attLib

def getCloudNetData(cloudNetFilePath,varName):
    
    cloudNetData = Dataset(cloudNetFilePath, 'r')
    height = cloudNetData['height'][:]
    time = cloudNetData['time'][:]
    variableData = cloudNetData[varName][:]
    print(time)
    print(height)
    print(variableData)
            
    return variableData, height, time
    


def fillGaps(timeRef, originalDF, colName):
    
    # fill the absence of values    
    timeData = originalDF.index
    filledDF = pd.DataFrame(index=timeRef, columns=[colName],
                            data=np.zeros(len(timeRef),np.int)-999)

    for timeId, timeStamp in enumerate(timeData):
    
        try:
            data = originalDF[0][timeStamp]                          
            filledDF[colName][timeStamp:timaData[timeId+1]] = data
    
        except:
            data = originalDF[0][timeStamp]
            filledDF[colName][timeStamp:]=data
    
        #print timeStamp, timeData[timeId+1]
        
    filledDF[filledDF<0]=0
    
    return filledDF


def getRainFlag(cloudNetFilePath, timeRef, colName,
                year, month, day):
    
    cloudNetCateg, heightData, timeData = getCloudNetData(cloudNetFilePath, 'category_bits')

    cloudNetCateg[cloudNetCateg>3]=-1
    cloudNetCateg[cloudNetCateg>0]=1
    cloudNetCateg[cloudNetCateg<0]=0

    timeData = attLib.convCloudTime2HumTime(timeData)

    timeDataCorr = pd.datetime(int(year), int(month), int(day))\
                               + pd.to_timedelta(np.array(timeData),
                               unit='s')

    categDF = pd.DataFrame(index=timeDataCorr,
                           columns=heightData,
                           data=cloudNetCateg)

    flagRain = pd.DataFrame(categDF.sum(axis=1))
    flagRain[flagRain>0]=1

    finalFlagRain = fillGaps(timeRef, flagRain, 'rainFlag')
    
    return finalFlagRain
    

def getHatProData(hatproFilePath,varName):
    
    cloudNetData = Dataset(hatproFilePath, 'r')
    
    units = getattr(cloudNetData.variables['time'],'units').split(' ')
    epoch = pd.to_datetime(' '.join([units[-3], units[-2]]))
    
    time = cloudNetData['time'][:]
    variableData = cloudNetData[varName][:]
    humanTime = epoch + pd.to_timedelta(time, unit='S' )
    
    return variableData, humanTime


def getLwpFlag(lwpFilePath, timeRef, colName,
               year, month, day):

    hatproLwp, timeData = getHatProData(lwpFilePath, 'clwvi')
    lwpDF = pd.DataFrame(index=timeData, data=hatproLwp*1000.)
    
    lwpDF['times']=timeData
    lwpDF=lwpDF.drop_duplicates(subset=['times'])
    del lwpDF['times']
    
    lwpDF = lwpDF.dropna(how='any')
    lwpDF[lwpDF <= 200] = 0
    lwpDF[lwpDF > 200] = 1

    finalFlagLwp= fillGaps(timeRef, lwpDF, colName)

    return finalFlagLwp

#def getLwpFlag(cloudNetFilePath, timeRef, colName, 
               #year, month, day):
    
#    cloudNetLwp, heightData, timeData = getCloudNetData(cloudNetFilePath, 'lwp')
#    timeData = attLib.convCloudTime2HumTime(timeData)
#    timeDataCorr = pd.datetime(int(year), int(month), int(day))\
               #                + pd.to_timedelta(np.array(timeData),
               #                unit='s')

#    lwpDF = pd.DataFrame(index=timeDataCorr, data=cloudNetLwp)
#    lwpDF = lwpDF.dropna(how='any')
#    lwpDF[lwpDF <= 200] = 0
#    lwpDF[lwpDF > 200] = 1
    
#    finalFlagLwp= fillGaps(timeRef, lwpDF, 'lwpFlag')
    
#    return finalFlagLwp


def getFlag(dataFrame, criteria): 

   dataFrame[np.isnan(dataFrame)==True]=0
   dataFrameFlag = dataFrame.copy()
   dataFrameAna = dataFrame.copy()
   dataFrameFlag[dataFrameAna >= criteria]=0
   dataFrameFlag[dataFrameAna < criteria]=1

   return dataFrameFlag


def getVarianceFlagOld(dt1, dt2):
    
    aveWind = 15*2
    variance =  (dt1 - dt2).rolling(aveWind, 
                                    min_periods=aveWind-2,
				    center=True).var()
    
    varianceFlag = variance.copy()
    varianceFlag[variance > 2] =1
    varianceFlag[variance <= 2] =0
    varianceFlag = varianceFlag.fillna(1)

    varianceFlagFinal = np.array(varianceFlag.values, np.uint16) << 13
    varianceFlagDF = pd.DataFrame(index=variance.index, columns=variance.columns,
                                  data=varianceFlagFinal)

    return(varianceFlagDF)


def getVarianceFlag(dt1, dt2):

    tempDS = xr.Dataset({'dt1':(['time', 'range'], dt1.values),
                         'dt2':(['time', 'range'], dt2.values)},
                         coords={'time':dt1.index,
                         'range':np.array(dt1.columns.values,float)}
                       )

    dwr = tempDS['dt1']-tempDS['dt2']
    var = dwr.rolling(time=30, min_periods=15, center=True).var()
    var.values = var.fillna(3)

    varFlag = var.copy()
    varFlag.values[var.values>2] = 1
    varFlag.values[var.values<=2] = 0
    varFlag.values = varFlag.fillna(1)

    varFlagFinal = np.array(varFlag.values, np.uint16) << 13
    varFlagDF = pd.DataFrame(index=dt1.index, columns=dt1.columns,
                             data=varFlagFinal)
    
    return(varFlagDF)




def getUnifiedFlag(rainFlagDF, lwpFlagDF,
                   corrFlagDF, pointFlagDF):

   rainFlag = rainFlagDF.loc[:,rainFlagDF.columns==108]
   
   rainFlag = np.array(rainFlag, np.uint16)<<7
   
   lwpFlag = lwpFlagDF.loc[:,lwpFlagDF.columns==108]
   lwpFlag = np.array(lwpFlag, np.uint16)<<6

   corrFlag = corrFlagDF.loc[:,corrFlagDF.columns==108]
   corrFlag = np.array(corrFlag, np.uint16)<<14

   pointFlag = pointFlagDF.loc[:,pointFlagDF.columns==108]
   pointFlag = np.array(pointFlag, np.uint16)<<15

   finalFlag = rainFlag + lwpFlag + corrFlag + pointFlag
   flagDF = pd.DataFrame(index=rainFlagDF.index, columns=['flag'],
                         data=finalFlag)

   return flagDF

