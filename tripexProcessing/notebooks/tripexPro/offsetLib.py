import numpy as np
from netCDF4 import Dataset
import pandas as pd
import tripexLib as trLib
import glob
import pandas as pd
import matplotlib.pyplot as plt


def getDataFrameList(fileList, variableDic):

    dataFrameList = []
    epoch = None	    

    varNames = variableDic.keys()
    for i, variable in enumerate(varNames):
        dataFrameList.append(pd.DataFrame())
        #print variable

    for fileName in fileList:
    
	#print fileName 
        rootgrpRe = Dataset(fileName, 'r')
        times = rootgrpRe.variables['time'][:]
        ranges = rootgrpRe.variables['range'][:]
        epoch = pd.to_datetime('1970-01-01 00:00:00')
        #epoch = trLib.getEpochTime(rootgrpRe, 'X')
        times = epoch+pd.to_timedelta(times, unit='s')
        
        for i, variable in enumerate(varNames):
    
            try:
                data = rootgrpRe.variables[variable][:]
                data = data + variableDic[variable]['offset']               
        
            except:
                if variable == 'halt':
                    data = np.ones((len(times), len(ranges)))*ranges
                else:
                    data = np.ones((len(times), len(ranges)))*np.nan
       		
		        
            dataFrame = pd.DataFrame(data=data, columns=ranges, index=times)
            
            dataFrameList[i] = dataFrameList[i].append(dataFrame)
        
        rootgrpRe.close()

    for i, variable in enumerate(varNames):
        dataFrameList[i]['times']=dataFrameList[i].index
        dataFrameList[i]=dataFrameList[i].drop_duplicates(subset='times')
        del dataFrameList[i]['times']
    
    return dataFrameList, epoch


def getMaskedDF(dataFrameList, variable,
                ZeKaMax, ZeKaMin, heightMin,
                offsetPair):

    varNames = variable.keys()
    kaIndex = varNames.index('Ka_DBZ_H')
    kaThresholdMax = dataFrameList[kaIndex]>=ZeKaMax
    kaThresholdMin = dataFrameList[kaIndex]<=ZeKaMin
    
    cloudKa = dataFrameList[kaIndex]/dataFrameList[kaIndex]

    for varNameCorr in offsetPair:
        
        varIndex = varNames.index(varNameCorr)
	
        dataFrameList[varIndex][cloudKa!=1]=np.nan
#        dataFrameList[varIndex].loc[:,dataFrameList[varIndex].columns 
#                                    < heightMin]=np.nan
        dataFrameList[varIndex][kaThresholdMax]=np.nan
        dataFrameList[varIndex][kaThresholdMin]=np.nan
        
    var0 = varNames.index(offsetPair[0])
    var1 = varNames.index(offsetPair[1])

    return dataFrameList[var0], dataFrameList[var1]

#import mkPlots as mkP
def getOffset(dataFrame, dataFrameRef,
              timesBegin, timesEnd,
	      radRef, rad):
    
    offSetAr = np.zeros(1440) + np.nan
    stdDiffAr = np.zeros(1440) + np.nan
    validPointsAr = np.zeros(1440) + np.nan
    stdZe_RefAr = np.zeros(1440) + np.nan
    correlAr = np.zeros(1440) + np.nan
    percentAr = np.zeros(1440) + np.nan

    for i in range(len(timesBegin)):
        newStart = timesBegin[i]
        newEnd = timesEnd[i]

        Ze_Rad = dataFrame[newStart:newEnd]
        Ze_Ref = dataFrameRef[newStart:newEnd]
        #print(np.isfinite(Ze_Rad.any()))
        Ze_Rad = np.ma.masked_invalid(Ze_Rad)
        Ze_Ref = np.ma.masked_invalid(Ze_Ref)    
        Ze_Diff = Ze_Rad - Ze_Ref	
        tempColDF = pd.DataFrame()
        #tempColDF['var1']=Ze_Ref.flatten()
        #tempColDF['var2']=Ze_Rad.flatten()
        var1Cor = Ze_Ref.flatten()#
        var2Cor = Ze_Rad.flatten()#
        var1Cor[~np.isfinite(var2Cor)] = np.nan
        var2Cor[~np.isfinite(var1Cor)] = np.nan
        var1Cor = var1Cor[np.isfinite(var1Cor)]
        var2Cor = var2Cor[np.isfinite(var2Cor)]
        tempColDF['var1']=var1Cor
        tempColDF['var2']=var2Cor
        correl = tempColDF['var1'].corr(tempColDF['var2'])
#       mkP.plotScatOff(Ze_Ref.flatten(), Ze_Rad.flatten(), 
#                      newStart, newEnd, radRef, rad, correl)###
        offset = np.mean(Ze_Diff)
        stdDiff = np.std(Ze_Diff)
        stdZe_Ref = np.std(Ze_Ref)   
        #correl = np.correlate(Ze_Ref.flatten(),Ze_Rad.flatten())
        if offset is np.ma.masked:
            offset = np.nan
        if stdDiff is np.ma.masked:
            stdDiff = np.nan    
        if stdZe_Ref is np.ma.masked:
            stdZe_Ref = np.nan
 
        invalidPoints = Ze_Diff[np.isnan(Ze_Diff)].shape[0]
        totalPoints = Ze_Diff.shape[0]*Ze_Diff.shape[1]
        validPoints = totalPoints - invalidPoints
        percent = (validPoints*100.) / float(totalPoints)
   
        if validPoints >= 300:	
            offSetAr[i] = offset
            correlAr[i] = correl
        else:
            offSetAr[i] = np.nan
            correlAr[i] = np.nan

        stdDiffAr[i] = stdDiff
        validPointsAr[i] = validPoints
        stdZe_RefAr[i] = stdZe_Ref
        percentAr[i] = percent


    return (offSetAr, stdDiffAr, validPointsAr,
            stdZe_RefAr, correlAr,percentAr)
        

def getParameterTimeSerie(parameters, timeFreq):
    
    repeatParameter = 60. / float(timeFreq[:-1]) 
    offsetTimeSerie = np.repeat(parameters[0],repeatParameter)
    stdDiffTimeSerie = np.repeat(parameters[1],repeatParameter)
    validPointsTimeSerie = np.repeat(parameters[2],repeatParameter)
    stdKaTimeSerie = np.repeat(parameters[3],repeatParameter)
    correlTimeSerie = np.repeat(parameters[4],repeatParameter)
    percentTimeSerie = np.repeat(parameters[5],repeatParameter)
    
    return (offsetTimeSerie, stdDiffTimeSerie, 
           validPointsTimeSerie, stdKaTimeSerie,
	   correlTimeSerie, percentTimeSerie)


def applyOffsetCorr(dataFrameToCorret, offset,
                    validPoints, thresholdPoints):
    
    offset[validPoints<thresholdPoints]=0
    dataFrameCorrected = dataFrameToCorret + offset
    
    return dataFrameCorrected


def getParamDF(parameter, timeRef, rangeRef):
        
    parameter[np.isnan(parameter)]=np.nan
    parameterDF=pd.DataFrame(index=timeRef, columns=rangeRef,
                          data=np.tile(parameter,(len(rangeRef),1)).T)
    
    return parameterDF


def getShiftedTemp(interpTempDF, timeRef, rangeRef):
    
    interpTempDF[interpTempDF > 0] = 1 
    interpTempDF[interpTempDF < 0] = -1

    tempArr = np.array(interpTempDF).T
    baseArr = np.zeros((40, tempArr.shape[1]))+1
    shiftedTempArr = np.vstack((baseArr, tempArr))
    shiftedTempArr = np.delete(shiftedTempArr,
                               range(shiftedTempArr.shape[0]-40,
                                     shiftedTempArr.shape[0]),0)

    shiftedTempArr[tempArr>0]=2
    shiftedTempDF = pd.DataFrame(index=timeRef, 
                                 columns=rangeRef, 
                                 data=shiftedTempArr.T)
    
    return shiftedTempDF



def temperatureMask(shiftedTempMaskDF, dataFrameListToMask,
                   variableToMask, timeRef, rangeRef):

    dfList = []
    shiftedTempMaskArr = np.array(shiftedTempMaskDF)

    for variable in variableToMask:
        
        dataArr = np.array(dataFrameListToMask[variableToMask.index(variable)])
        maskedData = np.ma.masked_where(shiftedTempMaskArr > 0, dataArr)
        
        maskedDataDF = pd.DataFrame(index=timeRef, columns=rangeRef, data=maskedData)
        
        dfList+=[maskedDataDF]
        
    return dfList


def getRangeBnds(rangeRef, rangeTolerance):
#    range_bnds = np.ones((len(rangeRef),2))
    range_bnds = np.ones((len(rangeRef),4))
    range_bnds[:,0:1] = (rangeRef - rangeTolerance).reshape(len(rangeRef),1)
    range_bnds[:,1:2] = (rangeRef - rangeTolerance).reshape(len(rangeRef),1)#
    range_bnds[:,2:3] = (rangeRef + rangeTolerance).reshape(len(rangeRef),1)
    range_bnds[:,3:4] = (rangeRef + rangeTolerance).reshape(len(rangeRef),1)#

    return range_bnds


def getTimeBnds(timeRef, timeTolerance):
#    time_bnds = np.ones((len(timeRef),2))
    time_bnds = np.ones((len(timeRef),4))
    time_bnds[:,0:1] = np.array(timeRef - pd.to_timedelta(2,'s'), float).reshape(len(timeRef),1)
    time_bnds[:,1:2] = np.array(timeRef + pd.to_timedelta(2,'s'), float).reshape(len(timeRef),1)
    time_bnds[:,2:3] = np.array(timeRef + pd.to_timedelta(2,'s'), float).reshape(len(timeRef),1)
    time_bnds[:,3:4] = np.array(timeRef - pd.to_timedelta(2,'s'), float).reshape(len(timeRef),1)

    time_bnds = time_bnds / 10**9
    
    return time_bnds



