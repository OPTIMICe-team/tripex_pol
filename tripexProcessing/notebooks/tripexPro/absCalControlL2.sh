#-----------------------------
#Shell script to control the
#resample process for entire 
#data set 
#
#--File Definition------------
inputCloudNetPath='/data/data_hatpro/jue/cloudnet/juelich/processed/categorize'
prefixCloudNet='juelich_categorize.nc'

inputPathL1='/home/jdias/Projects/radarDataResampled/absCalCorr'
prefixL1='tripex_3fr_L1_mom'

outputPathL2='/home/jdias/Projects/radarDataResampled/absCalCorr'
prefixL2='tripex_3fr_L2_mom'
#-----------------------------

#--Time Definition------------
#years=('2015' '2016')
years=('2016') #test
#---year
for year in ${years[@]}
do

if [ $year == '2015' ]
then
months=('08' '12')
#months=('11') #test
fi

if [ $year == '2016' ]
then
months=('08')
fi

#---month
for month in ${months[@]}
do

if [ $month == '08' ]
then
#days=('11' '12' '13' '14' '15' '16' '17' '18' '19' '20'\
#      '21' '22' '23' '24' '25' '26' '27' '28' '29' '30')
days=('11') #test
fi

if [ $month == '12' ]
then
#days=('01' '02' '03' '04' '05' '06' '07' '08' '09' '10'\
#      '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'\
#      '21' '22' '23' '24' '25' '26' '27' '28' '29' '30'\
#      '31')
days=('17')
fi

if [ $month == '01' ]
then
days=('01')
fi

#---day
for day in ${days[@]}
do

echo Date: $(tput setaf 1) $year $month $day $(tput sgr 0)
#---time


timeFreq=4s
timeTolerance=2s
#-----------------------------

#--Range Definition
beguinRangeRef=100
endRangeRef=12000
rangeFreq=30
rangeTolerance=17
#-----------------------------

#--Radar Offset Definitions---
#
zeKaMax=0
zeKaMin=-11
heightThreshold=5000
timeWindowLenght=5
pointsThreshold=300
#-----------------------------

#--Execution------------------
python tripexL2.py $inputCloudNetPath $prefixCloudNet\
                   $inputPathL1 $prefixL1 $outputPathL2 $prefixL2\
                   $year $month $day $timeFreq $timeTolerance \
                   $beguinRangeRef $endRangeRef $rangeFreq $rangeTolerance\
                   $zeKaMax $zeKaMin $heightThreshold $timeWindowLenght\
                   $pointsThreshold

#-----------------------------

done
done
done
