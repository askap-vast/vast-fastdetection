#!/bin/bash 

sbid=$1

path=`pwd`/SB$sbid
loc='/home/ymwang/vast_fastdetection' # code location 


SCRIPTS=$path/scripts


module use /home/app/modulefiles/
module load python/cpu-3.7.4


python $loc/tools/get_everything_ready.py $sbid $path

bash $SCRIPTS/download_selavy.sh
bash $SCRIPTS/bash_GETDATA.sh
bash $SCRIPTS/bash_CHECKDATA.sh

for i in {00..35}
#for i in {33..34}
do
    bash $loc/tools/submit_onebeam.sh beam$i $path &
    sleep 1s
done






