import glob
import numpy as np
import pandas as pd
from sys import argv

import xarray as xr
import writeData
import externalData as extLib
import offsetLib as offLib
import attenuationLib as attLib
import filters as filt
import qualityFlag as qFlag
import time as time1
import matplotlib.pyplot as plt
#--File Paths Definition------
#cloudNetFiles input
#dayIndices = argv[1]
#print(dayIndices.split(','))
#dayIndices = dayIndices.split(',')
#for dayIndex in dayIndices:
timeStart = time1.time()
#dayIndex = int(dayIndex)
    
    #-- define dateList
indStartEnd = argv[1]
indStart = indStartEnd.split(',')[0]
indEnd = indStartEnd.split(',')[1]
indStart = int(indStart)
indEnd = int(indEnd)
dateStart = pd.to_datetime('20190122'); dateEnd = pd.to_datetime('20190221')    
dateList = pd.date_range(dateStart, dateEnd,freq='D')

for dayIndex,date2proc in enumerate(dateList[indStart:indEnd]):
    #date2proc = dateList[dayIndex]; 
    date2end = dateList[dayIndex+indStart+1]
    print(date2proc)
    cloudNetPath = argv[2]#'/data/data_hatpro/jue/cloudnet/juelich/processed/categorize'#argv[1]
    cloudNetFileID =argv[3]#'juelich_categorize.nc'# argv[2]
    #tripexL1FilesL1 input
    dataPathRe = argv[4]#'/work/lvonterz/Cheops/output'#argv[3]
    fileId = argv[5]#'moments'#argv[4]
    #tripexL1FilesL2 output
    outputPath = argv[6]#'/work/lvonterz/tripexPro/outputtest'#argv[5]
    prefixL2 = argv[7]#'tripex_3fr_L2_mom'#argv[6]
    #-----------------------------

    #--Time Definitions-----------
    #year = argv[7]
    #month = argv[8]
    #day = argv[9]
    timeFreq= argv[8]
    timeTolerance = argv[9]

    #Define the reference time
    start = date2proc
    end = date2end - pd.offsets.Second(1)
    year = date2proc.strftime('%Y')
    month = date2proc.strftime('%m')
    day = date2proc.strftime('%d')
    print(year,month,day)
    timeRef = pd.date_range(start, end, freq=timeFreq)
    print(start,end,timeFreq)
    print(timeRef)
    #timesRef use the index from tripex data
    #-----------------------------

    #--Range Definitions----------
    beguinRangeRef = int(argv[10]) #botton
    endRangeRef = int(argv[11]) #top #original 12000
    rangeFreq = int(argv[12]) #rangeFreq
    rangeTolerance = int(argv[13]) #tol
    rangeRef = np.arange(beguinRangeRef, endRangeRef, rangeFreq)
    #rangeRef use the columns from tripex data
    #-----------------------------

    #--Radar Variables------------
    #Definitions to apply the offset correction
    #Ka_DBZ_HMax = int(argv[16]) #[dBZ]
    #Ka_DBZ_HMin = int(argv[17]) #[dBZ]
    heightThreshold = int(argv[14])#5000 #[m]
    timeWindowLenght = int(argv[15])#60 #[min]
    thresholdPoints = int(argv[16])#300 #[Threshold of Points]
    zeOffsetKa = float(argv[17])
    zeOffsetW = float(argv[18]) 
    zeOffsetX = float(argv[19])
    #------------------------------

    hatproPath = argv[20]#'/data/obs/site/jue/tophat/l2'
    hatproFileID =argv[21]#'sups_joy_mwr00_l2_clwvi_p00'
    
    #--Files to work--------------- 
    cloudNetFile = ('_').join([year+month+day,cloudNetFileID])
    cloudNetFilePath = ('/').join([cloudNetPath, year, cloudNetFile])

    hatproFile = ('_').join([hatproFileID, year+month+day+'*.nc'])
    hatproFilePath = ('/').join([hatproPath, year,month,day, hatproFile])
    #quit()
    hatproFileName = glob.glob(hatproFilePath)[0]

    #print cloudNetFilePath

    fileDate = ('').join([year, month, day])
    print(('/').join([dataPathRe,fileDate+'_'+fileId]))
    fileList = glob.glob(('/').join([dataPathRe,fileDate+'_'+fileId]))
    fileList = sorted(fileList)
    print(fileList)
    #-----------------------------

    #--Radar definitions----------
    #Radar Frequency
    radarFreqs = [9.4, 35.5, 95]#[GHz]

    #variables to be corrected
    variable={'X_DBZ_H':{'offset': zeOffsetX, 'colRange':(-25, 35), 'freq':9.4},
            'Ka_DBZ_H':{'offset': zeOffsetKa, 'colRange':(-25, 35), 'freq':35.5},
            'W_DBZ_H':{'offset': zeOffsetW, 'colRange':(-25, 35), 'freq':95}}
    varNames = variable.keys()
    #-----------------------------

    #--Attenuation correction-----
    results, time, height_M, temp, relHum, press = \
    attLib.getAtmAttPantra(cloudNetFilePath, radarFreqs)
    #print(results)
    interpAttDataList, qualityFlagList = \
    attLib.getInterpQualFlagList(results, time, timeRef, 
                                    timeTolerance, height_M,
                                    rangeRef, rangeTolerance, 
                                    radarFreqs, year, month, day)

    interpAttDataList = attLib.changeAttListOrder(interpAttDataList,
                                                variable, radarFreqs)
    qualityFlagList = attLib.changeAttListOrder(qualityFlagList, 
                                                variable, radarFreqs)
    #-----------------------------


    #--Copy temp from CLOUDNET----
    tempCel = temp - 273.15

    #print 'temp'
    #print tempCel[0]
    #print 'press'
    #print press[0]
    #print 'hum'
    #print relHum[0]


    resampledTemp = attLib.getResampledTimeRange(rangeRef, rangeTolerance, timeRef,
                                                time, timeTolerance, tempCel, year,
                                                month, day, height_M)

    interpTemp, qualityFlagTemp = attLib.getInterpData(time, timeRef, height_M,
                                                    resampledTemp, tempCel,
                                                    rangeRef)

    interpTempDF = pd.DataFrame(index=timeRef, columns=rangeRef, 
                                data=interpTemp.T)

    #-----------------------------

    #print interpTempDF[100]

    #--Copy press from CLOUDNET----
    resampledPress = attLib.getResampledTimeRange(rangeRef, rangeTolerance, timeRef,
                                                time, timeTolerance, press, year,
                                                month, day, height_M)

    interpPress, qualityFlagPress = attLib.getInterpData(time, timeRef, height_M,
                                                        resampledPress, press,
                                                        rangeRef)

    interpPressDF = pd.DataFrame(index=timeRef, columns=rangeRef, 
                                data=interpPress.T)

    #-----------------------------

    #print interpPressDF[100]

    #--Copy relHum from CLOUDNET----
    resampledRelHum = attLib.getResampledTimeRange(rangeRef, rangeTolerance, timeRef,
                                                time, timeTolerance, relHum, year,
                                                month, day, height_M)

    interpRelHum, qualityFlagRelHum = attLib.getInterpData(time, timeRef, height_M,
                                                        resampledRelHum, relHum,
                                                        rangeRef)

    interpRelHumDF = pd.DataFrame(index=timeRef, columns=rangeRef, 
                                data=interpRelHum.T)

    #-----------------------------
    #print interpRelHumDF[100]


    #--Offset correction----------here he adds the predefined offset!
    dataFrameList, epoch = offLib.getDataFrameList(fileList, variable)
    #dataFrameListNoCorrec, epoch = offLib.getDataFrameList(fileList, variable)


    #it removes extreme values from reflectively velocity 
    dataFrameList = filt.removeOutliersZeKa(dataFrameList, variable)
    #it removes the clutter from X band
    dataFrameList = filt.removeClutter(dataFrameList, variable, 'X_DBZ_H', 400)
    #it removes the clutter from Ka band
    dataFrameList = filt.removeClutter(dataFrameList, variable, 'Ka_DBZ_H', 400)
    #it removes the clutter from W band
    dataFrameList = filt.removeClutter(dataFrameList, variable, 'W_DBZ_H', 400)


    shiftedTempDF = interpTempDF.copy()
    shiftedTempDF = offLib.getShiftedTemp(shiftedTempDF, timeRef, rangeRef)

    #print(rangeRef)

    #Attenuation correction
    dataFrameListAtt = attLib.applyAttCorr(dataFrameList*1, interpAttDataList, variable)
    #print(dataFrameListAtt[0].values)
#    dsX = xr.DataArray(dataFrameListAtt[0].values,dims=('time','range'),coords={'time':timeRef,'range':rangeRef})
#    dsKa = xr.DataArray(dataFrameListAtt[1].values,dims=('time','range'),coords={'time':timeRef,'range':rangeRef})
#    dsW = xr.DataArray(dataFrameListAtt[2].values,dims=('time','range'),coords={'time':timeRef,'range':rangeRef})
#    ds = xr.Dataset({'X':dsX,'Ka':dsKa,'W':dsW})
#    ds.to_netcdf('/work/lvonterz/tripexProcessing/'+year+month+day+'_att_corr_no_off.nc')
#    print(ds)
#    '''
    #offset 
    dataFrameListMasked = attLib.applyAttCorr(dataFrameList*1, interpAttDataList, 
                                            variable)
    
    timeWindow = pd.to_timedelta( timeWindowLenght, unit='m')
    timesBegin = pd.date_range(start-timeWindow, end-timeWindow, freq='1min')
    timesEnd = pd.date_range(start+timeWindow, end+timeWindow, freq='1min')

    #offset X Ka
    offsetPairXKa = ['X_DBZ_H','Ka_DBZ_H']
    dataFrameListToXKa = attLib.applyAttCorr(dataFrameList*1, interpAttDataList, 
                                            variable)

    dataFrameListMaskedXKa = offLib.getMaskedDF(dataFrameListToXKa, variable, 
                                                0, -15, heightThreshold,
                                                offsetPairXKa)

    maskedTempDFlistXKa = offLib.temperatureMask(shiftedTempDF,
                                                dataFrameListMaskedXKa, 
                                                offsetPairXKa, timeRef,
                                                rangeRef)

    dataFrame = maskedTempDFlistXKa[offsetPairXKa.index('X_DBZ_H')]
    dataFrameRef = maskedTempDFlistXKa[offsetPairXKa.index('Ka_DBZ_H')]
    parametersXKa= offLib.getOffset(dataFrame, dataFrameRef,
                                    timesBegin, timesEnd, 
                                    'Ka', 'X')
    #print(parametersXKa[4][np.isfinite(parametersXKa[4])])
    #offset Ka W
    offsetPairKaW = ['Ka_DBZ_H','W_DBZ_H']
    dataFrameListToKaW = attLib.applyAttCorr(dataFrameList*1, interpAttDataList, 
                                            variable)

    dataFrameListMaskedKaW = offLib.getMaskedDF(dataFrameListToKaW, variable, 
                                                -10, -30, heightThreshold,
                                                offsetPairKaW)

    maskedTempDFlistKaW = offLib.temperatureMask(shiftedTempDF,
                                                dataFrameListMaskedKaW, 
                                                offsetPairKaW, timeRef,
                                                rangeRef)
    
    dataFrame = maskedTempDFlistKaW[offsetPairKaW.index('W_DBZ_H')]
    dataFrameRef = maskedTempDFlistKaW[offsetPairKaW.index('Ka_DBZ_H')]
    parametersWKa= offLib.getOffset(dataFrame, dataFrameRef,
                                    timesBegin, timesEnd,
                                    'Ka', 'W')
    parametersXKaTS = offLib.getParameterTimeSerie(parametersXKa, timeFreq)
    parametersWKaTS = offLib.getParameterTimeSerie(parametersWKa, timeFreq)

    ###(I mutipled the offset by -1 ) ( ## testing the offset)
    percentDfWKa = offLib.getParamDF(parametersWKaTS[5]*1, timeRef, rangeRef)
    offsetWKaDF = offLib.getParamDF(parametersWKaTS[0]*1, timeRef, rangeRef)
    validPointWKaDF = offLib.getParamDF(parametersWKaTS[2]*1, timeRef, rangeRef)
    correlXKaDF = offLib.getParamDF(parametersXKaTS[4]*1, timeRef, rangeRef)
    
    offsetWKaDF = offsetWKaDF*(-1)
    dataFrameListAtt[varNames.index('W_DBZ_H')] = \
    offLib.applyOffsetCorr(dataFrameListAtt[varNames.index('W_DBZ_H')],
                            offsetWKaDF*1, validPointWKaDF, 
                            thresholdPoints)

    percentDfXKa = offLib.getParamDF(parametersXKaTS[5]*1, timeRef, rangeRef)
    offsetXKaDF = offLib.getParamDF(parametersXKaTS[0]*1, timeRef, rangeRef)
    validPointXKaDF = offLib.getParamDF(parametersXKaTS[2]*1, timeRef, rangeRef)
    correlWKaDF = offLib.getParamDF(parametersWKaTS[4]*1, timeRef, rangeRef)

    offsetXKaDF = offsetXKaDF*(-1)
    dataFrameListAtt[varNames.index('X_DBZ_H')] = \
    offLib.applyOffsetCorr(dataFrameListAtt[varNames.index('X_DBZ_H')],
                            offsetXKaDF*1, validPointXKaDF, 
                            thresholdPoints)



    #quit()

    #--Quality Flags--------------

    rainFlag = qFlag.getRainFlag(cloudNetFilePath, timeRef, 'rainFlag',
                                year, month, day)
    rainFlagDF = offLib.getParamDF(rainFlag['rainFlag'], timeRef, rangeRef)
    lwpFlag = qFlag.getLwpFlag(hatproFileName, timeRef, 'lwpFlag',
                            year, month, day)
    lwpFlagDF = offLib.getParamDF(lwpFlag['lwpFlag'], timeRef, rangeRef)


    #offset flag X band
    valPoinFlagXKaDF = qFlag.getFlag(validPointXKaDF, 300)
    corrFlagXKaDF = qFlag.getFlag(correlXKaDF, 0.70)

    varFlagXKaDF = qFlag.getVarianceFlag(dataFrameListAtt[varNames.index('X_DBZ_H')],
                                        dataFrameListAtt[varNames.index('Ka_DBZ_H')])

    finalFlagXKa = qFlag.getUnifiedFlag(rainFlagDF, lwpFlagDF,
                                        corrFlagXKaDF, valPoinFlagXKaDF)

    finalFlagXKaDF = offLib.getParamDF(finalFlagXKa['flag'], timeRef, rangeRef) 
    finalFlagXKaDF = finalFlagXKaDF + varFlagXKaDF


    #offset flag W band
    valPoinFlagWKaDF = qFlag.getFlag(validPointWKaDF, 300)
    corrFlagWKaDF = qFlag.getFlag(correlWKaDF, 0.70)

    varFlagKaWDF = qFlag.getVarianceFlag(dataFrameListAtt[varNames.index('Ka_DBZ_H')],
                                        dataFrameListAtt[varNames.index('W_DBZ_H')])

    finalFlagWKa = qFlag.getUnifiedFlag(rainFlagDF, lwpFlagDF,
                                        corrFlagWKaDF, valPoinFlagWKaDF)

    finalFlagWKaDF = offLib.getParamDF(finalFlagWKa['flag'], timeRef, rangeRef)
    finalFlagWKaDF = finalFlagWKaDF + varFlagKaWDF
    print(finalFlagWKaDF.dtypes)
    #quit()
    #-----------------------------



    #--radar definitions --------

    coordenates = {'lat':{'data':50.9086},
                'lon':{'data':6.4135},
                'zsl':{'data':112.5}
    }



    externalData = {'freq_sb_x':{'data':np.array(9.4*10**9,np.float32)},
                    'freq_sb_ka':{'data':np.array(35.5*10**9,np.float32)},
                    'freq_sb_w':{'data':np.array(94*10**9,np.float32)},
                    'radar_beam_width_x':{'data':np.array(1.3,np.float32)},
                    'radar_beam_width_ka':{'data':np.array(0.6,np.float32)},
                    'radar_beam_width_w':{'data':np.array(0.5,np.float32)},
    }


    #pd.to_timedelta(2,)


    bnds = {'time_bnds':{'data': offLib.getTimeBnds(timeRef, timeTolerance)},
            'range_bnds':{'data':offLib.getRangeBnds(rangeRef, rangeTolerance)},
    }

    #-----------------------------


    #--Copy data from L1----------
    data = None
    variableToCopy={'X_VEL_H':{'data':data, 'offset':0, 'outName':'X_VEL_H'},
                    'Ka_VEL_H':{'data':data, 'offset':0, 'outName':'Ka_VEL_H'},
                    'W_VEL_H':{'data':data, 'offset':0, 'outName':'W_VEL_H'},
                    'X_WIDTH_H':{'data':data, 'offset':0, 'outName':'X_WIDTH_H'},
                    'Ka_WIDTH_H':{'data':data, 'offset':0, 'outName':'Ka_WIDTH_H'},
                    'W_WIDTH_H':{'data':data, 'offset':0, 'outName':'W_WIDTH_H'},
                    'X_SK_H':{'data':data, 'offset':0, 'outName':'X_SK_H'},
                    'Ka_SK_H':{'data':data, 'offset':0, 'outName':'Ka_SK_H'},
                    'W_SK_H':{'data':data, 'offset':0, 'outName':'W_SK_H'},                    
                    #'LDRg_Ka':{'data':data, 'offset':0, 'outName':'ldr_ka'},
                    #'KDP_X':{'data':data, 'offset':0, 'outName':'kdp_x'},
                    #'PhiDP_X':{'data':data, 'offset':0, 'outName':'phidp_x'},
                    #'RhoHV_X':{'data':data, 'offset':0, 'outName':'rhohv_x'},
                    #'ZDR_X':{'data':data, 'offset':0, 'outName':'zdr_x'},
                    }

    dataCopiedDFList, epoch = offLib.getDataFrameList(fileList, 
                                                    variableToCopy)
    #it removes the clutter from X band
    dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'X_VEL_H', 400)
    dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'X_WIDTH_H', 400)
    dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'X_SK_H', 400)
    #dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'KDP_X', 700)
    #dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'PhiDP_X', 700)
    #dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'RhoHV_X', 700)
    #dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'ZDR_X', 700)

    #it removes the clutter from Ka band
    dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'Ka_VEL_H', 400)
    dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'Ka_WIDTH_H', 400)
    dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'Ka_SK_H', 400)
    #dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'LDRg_Ka', 400)


    #it removes the clutter from W band
    dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'W_VEL_H', 400)
    dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'W_WIDTH_H', 400)
    dataFrameList = filt.removeClutter(dataCopiedDFList, variableToCopy, 'W_SK_H', 400)

    variableToCopyTemp={}
    for variable in variableToCopy.keys():

        outName = variableToCopy[variable]['outName']
        variableToCopyTemp[outName]={'data':data}

    for variable in variableToCopy.keys():
        
        varNamesToCopy = variableToCopy.keys()
        outName = variableToCopy[variable]['outName']
        variableToCopyTemp[outName]['data'] = \
            dataCopiedDFList[varNamesToCopy.index(variable)]

    variableToCopy = variableToCopyTemp

    #-----------------------------

    #from IPython.core.debugger import Tracer ; Tracer()()
    #--Write data-----------------

    variableOutPut={'X_DBZ_H':{'data':dataFrameListAtt[varNames.index('X_DBZ_H')]},
                    'Ka_DBZ_H':{'data':dataFrameListAtt[varNames.index('Ka_DBZ_H')]},
                    'W_DBZ_H':{'data':dataFrameListAtt[varNames.index('W_DBZ_H')]},
                    'pia_x':{'data':interpAttDataList[varNames.index('X_DBZ_H')]},
                    'pia_ka':{'data':interpAttDataList[varNames.index('Ka_DBZ_H')]},
                    'pia_w':{'data':interpAttDataList[varNames.index('W_DBZ_H')]},
                    'offset_x':{'data':offsetXKaDF},
                    'offset_w':{'data':offsetWKaDF},
                # 'valDat_x':{'data':validPointXKaDF},
                # 'valDat_w':{'data':validPointWKaDF},
                # 'correlation_X':{'data':correlXKaDF},
                # 'correlation_W':{'data':correlWKaDF},
                # 'corrFlag_x':{'data':corrFlagXKaDF},
                # 'pointFlag_x':{'data':valPoinFlagXKaDF},
                # 'corrFlag_w':{'data':corrFlagWKaDF},
                # 'pointFlag_w':{'data':valPoinFlagWKaDF},
                # 'rainFlag_x':{'data':rainFlagDF},
                # 'lwpFlag_x':{'data':lwpFlagDF},
                    'pa':{'data':interpPressDF},	
                    'hur':{'data':interpRelHumDF},
                    'ta':{'data':interpTempDF},
                    'quality_flag_offset_x':{'data':finalFlagXKaDF},
                    'quality_flag_offset_w':{'data':finalFlagWKaDF},
                }

    #for indexWrite, timeStart in enumerate(timesBeginWrite):
            
    dateName = start.strftime('%Y%m%d')
    #print(dateName)    
    outPutFile = ('_').join([dateName,prefixL2+'.nc'])
    outPutFilePath = ('/').join([outputPath, outPutFile])
    #print(outPutFilePath)
    timeRefUnixWrt = np.array(timeRef, float)
    timeRefUnixWrt = timeRefUnixWrt/10.**9
    
    rootgrpOut = writeData.createNetCdf(outPutFilePath, prefixL2)


    for varNameOut in sorted(coordenates.keys()):
        
        varListName = varNameOut.split('_')
        if len(varListName) > 1:
            varFinalName = '_'.join(varListName[:-1])
            sensor = varListName[-1]
        else:
            varFinalName = varListName[0]
            sensor = ''        
        dataDF = coordenates[varNameOut]['data']
        dataToWrite = np.array(dataDF)
        var_Written = writeData.createOneValvariable(rootgrpOut, dataToWrite,
                                                    varFinalName, sensor, prefixL2)

    #for varNameOut in sorted(bnds.keys()):
            
    #    dimName, varName = varNameOut.split('_')
    #    dataToWrite = bnds[varNameOut]['data']
    #    var_Written = writeData.createBndsVariable(rootgrpOut, dataToWrite,
    #                                               varNameOut, dimName)
    
    nv_dim = writeData.createNvDimension(rootgrpOut, prefixL2)
    
    time_ref = writeData.createTimeDimension(rootgrpOut, timeRefUnixWrt, prefixL2)
    dataToWrite = bnds['time_bnds']['data']
    var_Written = writeData.createBndsVariable(rootgrpOut, dataToWrite,
                                            'time_bnds', 'time')
    
    range_ref = writeData.createRangeDimension(rootgrpOut, rangeRef, prefixL2)
    dataToWrite = bnds['range_bnds']['data']
    var_Written = writeData.createBndsVariable(rootgrpOut, dataToWrite,
                                            'range_bnds', 'range')


    for varNameOut in sorted(variableOutPut.keys()):
        if (varNameOut == 'X_DBZ_H') or (varNameOut == 'W_DBZ_H'): 
            radar = varNameOut[0]
            varFinalName = 'DBZ' 
        elif varNameOut == 'Ka_DBZ_H':
            radar = 'Ka'
            varFinalName = 'DBZ'
        else:
            varListName = varNameOut.split('_')
        
            if len(varListName) > 1:
                varFinalName = '_'.join(varListName[:-1])
                radar = varListName[-1]

            else:
                varFinalName = varListName[0]
                radar = ''        
        print(varFinalName)
        print(radar)     
        dataDF = variableOutPut[varNameOut]['data']
   #     varFinalName = varNameOut # das muss dann weg wenn wir bindestrich haben
        if varFinalName == 'quality_flag_offset':
            dataToWrite = np.array(dataDF.astype(np.uint16))
            var_Written = writeData.createVariable(rootgrpOut, dataToWrite,
                                                varFinalName, varNameOut,
                                                radar, prefixL2, np.uint16)
        else:
            dataToWrite = np.array(dataDF.astype(np.float32))
            var_Written = writeData.createVariable(rootgrpOut, dataToWrite,
                                                varFinalName, varNameOut,
                                                radar, prefixL2, np.float32)
    #It writes the data from L1 in L2 file
    for varNameOut in sorted(variableToCopy.keys()):

        #it removes the noise from v_Ka
        if varNameOut == 'bla': #'Ka_VEL_H': # we already removed the noise!!
                
            #indexV_Ka = variableToCopy.keys().index('Ka_VEL_H')
            #indexKa_DBZ_H = varNames.index('Ka_DBZ_V')    
            variableToCopy['Ka_VEL_H']['data'] = \
                    filt.removeVelNoiseKa(variableOutPut['Ka_DBZ_H']['data'],
                                        variableToCopy['Ka_VEL_H']['data'])

        elif (varNameOut == 'kdp_x') or (varNameOut == 'phidp_x')\
            or (varNameOut == 'rhohv_x') or (varNameOut == 'zdr_x'):

            variableToCopy[varNameOut]['data'] = \
                    filt.removeVelNoiseKa(variableToCopy['X_VEL_H']['data'],
                                        variableToCopy[varNameOut]['data'])
                    #filt.removeVelNoiseKa(variableOutPut['X_VEL_H']['data'],
                    #                     variableToCopy[varNameOut]['data'])





        else:
            pass
        #-------------------------------       
        if (varNameOut=='X_VEL_H') or (varNameOut=='W_VEL_H') or (varNameOut=='Ka_VEL_H') or (varNameOut=='X_WIDTH_H') or (varNameOut=='W_WIDTH_H') or (varNameOut=='Ka_WIDTH_H') or (varNameOut=='Ka_SK_H') or (varNameOut=='X_SK_H') or (varNameOut=='W_SK_H'): 
            #radar = varNameOut[0]
            radar,varFinalName,orient = varNameOut.split('_')
            print(radar,varFinalName,orient)
            #varFinalName = varNameOut
        #elif (varNameOut == 'Ka_VEL_H') or (varNameOut == 'Ka_WIDTH_H'):
            #radar = 'Ka'
            #varFinalName = varNameOut
        else:
            print(varNameOut)
            varFinalName,radar = varNameOut.split('_')
        
        
        dataDF = variableToCopy[varNameOut]['data'] 
        dataToWrite = np.array(dataDF.astype(np.float32))
        var_Written = writeData.createVariable(rootgrpOut, dataToWrite,
                                                varFinalName, varNameOut,
                                                radar, prefixL2, np.float32)
            

    for varNameOut in sorted(externalData.keys()):

        varListName = varNameOut.split('_')
        if len(varListName) > 1:
            varFinalName = '_'.join(varListName[:-1])
            sensor = varListName[-1]
        else:
            varFinalName = varListName[0]
            sensor = ''        
        dataDF = externalData[varNameOut]['data']
        dataToWrite = np.array(dataDF)
        var_Written = writeData.createOneValvariable(rootgrpOut, dataToWrite,
                                                    varFinalName, sensor, prefixL2)

    

    rootgrpOut.close()
    timeEnd = time1.time()
    print(timeEnd-timeStart)

    #-----------------------------



    
