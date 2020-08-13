#-----------------------------
#Shell script to control the
#resample process for entire 
#data set 
#
#--File Definition------------
inputCloudNetPath='/data/data_hatpro/jue/cloudnet/juelich/processed/categorize'
prefixCloudNet='juelich_categorize.nc'

#inputPathL1='/home/jdias/Projects/radarDataResampled/data'
#inputPathL1='/data/optimice/tripex/tripex_level_01'
inputPathL1='/data/optimice/tripex/tripex_level_01_X_pol'

prefixL1='tripex_3fr_L1_mom'

#outputPathL2='/work2/tripex/samd'
#outputPathL2='/data/optimice/tripex/tripex_level_02_NOSENS'
#outputPathL2='/data/optimice/tripex/tripex_level_02_samd'
#outputPathL2='/data/optimice/tripex/tripex_level_02_test'
#outputPathL2='/data/optimice/tripex/testOffset/tempData'
#outputPathL2='/data/optimice/tripex/testOffset/newOffsetRegimeKa'
outputPathL2='/data/optimice/tripex/tripex_level_02_X_pol'
prefixL2='tripex_joy_tricr00_l2_any_v00'


hatproPath='/data/hatpro/jue/data/level2'
hatproFileID='sups_joy_mwr00_l2_clwvi_v01' 

#-----------------------------

#--Time Definition------------
years=('2015' '2016')
#years=('2015') #test
#---year
for year in ${years[@]}
do

if [ $year == '2015' ]
then
months=('11' '12')
#months=('12') #test
#months=('11') #test
fi

if [ $year == '2016' ]
then
months=('01')
fi

#---month
for month in ${months[@]}
do


if [ $month == '08' ]
then
days=('27') #test
fi


if [ $month == '11' ]
then
days=('11' '12' '13' '14' '15' '16' '17' '18' '19' '20'\
      '21' '22' '23' '24' '25' '26' '27' '28' '29' '30')
#days=('19') #test
#days=('24') #test
fi

if [ $month == '12' ]
then
days=('01' '02' '03' '04' '05' '06' '07' '08' '09' '10'\
      '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'\
      '21' '22' '23' '24' '25' '26' '27' '28' '29' '30'\
      '31')
#days=('04' '08' '21' '31')
#days=('07')
fi

if [ $month == '01' ]
then
days=('01' '02' '03' '04')
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
timeWindowLenght=15
pointsThreshold=300
zeOffsetKa=3.6 #4dB from T-matrix calculation
#-----------------------------

#--Execution------------------
python tripexL2Samd.py $inputCloudNetPath $prefixCloudNet\
                   $inputPathL1 $prefixL1 $outputPathL2 $prefixL2\
                   $year $month $day $timeFreq $timeTolerance \
                   $beguinRangeRef $endRangeRef $rangeFreq $rangeTolerance\
                   $zeKaMax $zeKaMin $heightThreshold $timeWindowLenght\
                   $pointsThreshold $zeOffsetKa $hatproPath $hatproFileID

#-----------------------------

done
done
done
