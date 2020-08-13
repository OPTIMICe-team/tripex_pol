#-----------------------------
#Shell script to control the
#resample process for a single file
#
#--File Definition------------
inputPath='/data/obs/site/jue/joyrad35/' #'/home/jdias/Projects/radarData'
outputPath='/home/lvonterz/tripexPro/outputtest/' # '/home/jdias/Pictures/2015112419_1729_3rd/data'
prefix='tripex_3fr_L1_mom'
#-----------------------------

#--Time Definition------------
year=2018
month=12
day=09
beguinTime=20
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
#radar='X'
#rangeOffSet=-15.5
#variables=('Ze' 'vd')
#variables=('Ze')
#for variable in ${variables[@]}
#do
#python tripexL1.py $inputPath $outputPath $prefix $year $month \
#                    $day $beguinTime $timeFreq $timeTolerance \
#                    $beguinRangeRef $endRangeRef $rangeFreq \
#                    $rangeTolerance $radar $rangeOffSet $variable
#done
#-----------------------------
#echo Radar: $(tput setaf 3) $radar Done $(tput sgr 0)


#--W band Radar setup---------
#
#radar='W'
#rangeOffSet=0
#variables=('Ze' 'vm' 'sigma')
#for variable in ${variables[@]}
#do
#python tripexL1.py $inputPath $outputPath $prefix $year $month \
#                    $day $beguinTime $timeFreq $timeTolerance \
#                    $beguinRangeRef $endRangeRef $rangeFreq \
#                    $rangeTolerance $radar $rangeOffSet $variable
#done
#-----------------------------
#echo Radar: $(tput setaf 3) $radar Done $(tput sgr 0)


#--Ka band Radar setup---------
#
radar='Ka'
rangeOffSet=0
variables=('Zg') # 'VELg' 'RMS' 'LDR')
for variable in ${variables[@]}
do
python tripexL1.py $inputPath $outputPath $prefix $year $month \
                      $day $beguinTime $timeFreq $timeTolerance \
                      $beguinRangeRef $endRangeRef $rangeFreq \
                      $rangeTolerance $radar $rangeOffSet $variable
done
#-----------------------------
echo Radar: $(tput setaf 3) $radar Done $(tput sgr 0)


