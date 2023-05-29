#!/bin/bash 

sbid=$1

path=`pwd`
loc='~/vast_fastdetection' # code location 
SCRIPTS=$path/SCRIPTS
LOGS=$path/LOGS

VIS=$path/SB$sbid/visibilities
MODELS=$path/SB$sbid/models
IMGS=$path/SB$sbid/images

mkdir $SCRIPTS
mkdir $LOGS
mkdir $path/SB$sbid
mkdir $VIS
mkdir $MODELS
mkdir $IMGS


module use /home/app/modulefiles/
module load python/cpu-3.7.4


python $loc/tools/get_everything_ready.py $sbid $SCRIPTS

cd $VIS
bash $SCRIPTS/download_selavy.sh


for i in {00..35}
do
    bash $SCRIPTS/submit_onebeam.sh beam$i &
    sleep 1s
done






