#!/bin/bash
years=('2015' '2016')

#---year
for year in ${years[@]}
do

if [ $year == '2015' ]
then
months=('11' '12')
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
      '21' '22' '23' '24' '25' '26' '27' '28' '29' '30'\
      '31')
fi

if ($month == '12')
then
days=('01' '02' '03' '04' '05' '06' '07' '08' '09' '10'\
      '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'\
      '21' '22' '23' '24' '25' '26' '27' '28' '29' '30'\
      '31')
fi

if [ $month == '12' ]
then
days=('01' '02' '03' '04')
fi

#---day
for day in ${days[@]}
do
beguinTimes=('00' '01' '02' '03' '04' '05' '06' '07' '08' '09' '10'\
             '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21'\
             '22' '23' )

echo Date: $(tput setaf 1) $year $month $day $(tput sgr 0)
#---time
for beguinTime in ${beguinTimes[@]}
do
echo $year, $month, $day, $beguinTime

done
done
done
done



#if [ $numero == '2015' ]
#then
#echo 'True'
#else
#echo 'False'
#fi
