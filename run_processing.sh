#!/bin/bash


#sbid=32043
sbid=$1
field=$2

module use /home/app/modulefiles/
module load python/cpu-3.7.4

python /home/ymwang/casda_test/casda_download.py $sbid visibility_download_"$sbid".sh


bash visibility_download_"$sbid".sh
sleep 5m

# check if all of visibility are downloaded without problem
bash visibility_download_"$sbid".sh
sleep 5m

# untar the downloaded data 
ls *"$sbid"*tar | xargs -n 1 -P 10 --replace tar xvf {}

echo "Download and untar data finished. "


module load python/cpu-3.6.5

cd /o9000/ASKAP/VAST/fast_test/vast-fastimager

python run_ALL_ym.py $field $sbid > logfiles/run_"$sbid".out 

echo "Submit all jobs. "



