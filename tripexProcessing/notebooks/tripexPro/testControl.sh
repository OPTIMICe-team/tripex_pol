#-----------------------------
#Shell script to control the
#resample process 
#
#--File Definition------------
inputPath='/home/jdias/Projects/radarData'
outputPath='/home/jdias/Projects/radarDataResampled'
prefix='tripex_3fr_L1_momTest'
#-----------------------------

#--Time Definition------------
year=2015
month=11
day=24
beguinTime=17
timeFreq=4s
timeTolerance=2s
#-----------------------------

#--Range Definition
beguinRangeRef=100
endRangeRef=12000
rangeFreq=30
rangeTolerance=17
#-----------------------------


#--W band Radar setup---------
#
radar='Ka'
rangeOffSet=2
variables=('Zg')
for variable in ${variables[@]}
do
python tripexL1.py $inputPath $outputPath $prefix $year $month \
                    $day $beguinTime $timeFreq $timeTolerance \
                    $beguinRangeRef $endRangeRef $rangeFreq \
                    $rangeTolerance $radar $rangeOffSet $variable
done
#-----------------------------
echo Radar: $(tput setaf 3) $radar Done $(tput sgr 0)



