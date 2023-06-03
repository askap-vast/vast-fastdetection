#!/bin/bash 

sbid=$1

path=`pwd`
loc='~/vast_fastdetection' # code location 


VIS=$path/SB$sbid/visibilities
MODELS=$path/SB$sbid/models
IMGS=$path/SB$sbid/images
CAND=$path/SB$sbid/candidates

SCRIPTS=$path/SB$sbid/SCRIPTS
LOGS=$path/SB$sbid/LOGS


mkdir $path/SB$sbid
mkdir $VIS
mkdir $MODELS
mkdir $IMGS
mkdir $CAND

mkdir $SCRIPTS
mkdir $LOGS


module use /home/app/modulefiles/
module load python/cpu-3.7.4


python $loc/tools/get_everything_ready.py $sbid $SCRIPTS

# cd $VIS
# bash $SCRIPTS/download_selavy.sh


# for i in {00..35}
# do
#     bash $SCRIPTS/submit_onebeam.sh beam$i &
#     sleep 1s
# done






