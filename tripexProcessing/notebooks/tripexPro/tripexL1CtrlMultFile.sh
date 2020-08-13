#-----------------------------
#Shell script to control the
#resample process for entire 
#data set 
#
#--File Definition------------
inputPath='/work/radarData'
#outputPath='/data/optimice/tripex/tripex_level_01'
outputPath='/data/optimice/tripex/tripex_level_01_X_pol'
#outputPath='/work2/tripex/tripex_level_01_X_pol'
prefix='tripex_3fr_L1_mom'
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
#months=('11') #test
fi

if [ $year == '2016' ]
then
months=('01')
fi

#---month
for month in ${months[@]}
do

if [ $month == '11' ]
then
days=('11' '12' '13' '14' '15' '16' '17' '18' '19' '20'\
      '21' '22' '23' '24' '25' '26' '27' '28' '29' '30')
#days=('24') #test
fi

if [ $month == '12' ]
then
days=('01' '02' '03' '04' '05' '06' '07' '08' '09' '10'\
      '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'\
      '21' '22' '23' '24' '25' '26' '27' '28' '29' '30'\
      '31')
#days=('04')
fi

if [ $month == '01' ]
then
days=('01' '02' '03' '04')
fi

#---day
for day in ${days[@]}
do
beguinTimes=('00' '01' '02' '03' '04' '05' '06' '07' '08' '09' '10'\
             '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21'\
             '22' '23' )
#beguinTimes=('00' '01' '02' '03' '04' '05' '06' '07' '08' '09' '10'\
#             '11' '12' )
#
#beguinTimes=('08') #test

echo Date: $(tput setaf 1) $year $month $day $(tput sgr 0)
#---time
for beguinTime in ${beguinTimes[@]}
do

timeFreq=4s
timeTolerance=2s
#-----------------------------

#--Range Definition
beguinRangeRef=100
endRangeRef=12000
rangeFreq=30
rangeTolerance=17
#-----------------------------

#--X band Radar setup---------
#
radar='X'
rangeOffSet=-15.4
variables=('Ze' 'vd' 'W' 'ZDR' 'KDP' 'PhiDP' 'RhoHV')
zeOffset=0
for variable in ${variables[@]}
do
python tripexL1.py $inputPath $outputPath $prefix $year $month \
                    $day $beguinTime $timeFreq $timeTolerance \
                    $beguinRangeRef $endRangeRef $rangeFreq \
                    $rangeTolerance $radar $rangeOffSet $variable\
		    $zeOffset

done
#-----------------------------
echo Radar: $(tput setaf 3) $radar Done $(tput sgr 0)


#--W band Radar setup---------
#
radar='W'
rangeOffSet=0
variables=('Ze' 'vm' 'sigma')
zeOffset=0
for variable in ${variables[@]}
do
python tripexL1.py $inputPath $outputPath $prefix $year $month \
                    $day $beguinTime $timeFreq $timeTolerance \
                    $beguinRangeRef $endRangeRef $rangeFreq \
                    $rangeTolerance $radar $rangeOffSet $variable\
		    $zeOffset
done
#-----------------------------
echo Radar: $(tput setaf 3) $radar Done $(tput sgr 0)


#--Ka band Radar setup---------
#
radar='Ka'
rangeOffSet=2.2
variables=('Ze' 'VELg' 'RMS' 'LDRg')
zeOffset=5 #dB (MITEK correction)
for variable in ${variables[@]}
do
python tripexL1.py $inputPath $outputPath $prefix $year $month \
                      $day $beguinTime $timeFreq $timeTolerance \
                      $beguinRangeRef $endRangeRef $rangeFreq \
                      $rangeTolerance $radar $rangeOffSet $variable\
		      $zeOffset
done
#-----------------------------
echo Radar: $(tput setaf 3) $radar Done $(tput sgr 0)

done
done
done
done
