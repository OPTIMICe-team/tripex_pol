#! /bin/bash -l

# send days to processing, using separate nodes for each day

#python create_day_list.py # create a list with days to run, index needs to be provided in order to know which day to use

inputdays=$1

sbatch --array=$inputdays batch_process_day.sh  




