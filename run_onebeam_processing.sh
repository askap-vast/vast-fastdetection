#!/bin/bash 

# bash run_onebeam_processing.sh 45965 24

sbid=$1
i=$2

path=`pwd`/SB$sbid
loc=$(dirname "$0") # code location 

echo "The script you are running has:"
echo "basename: [$(basename "$0")]"
echo "dirname : [$(dirname "$0")]"
echo "pwd     : [$(pwd)]"
echo
echo The output will store to $path
echo The scripts are located in $loc


SCRIPTS=$path/scripts
DATA=$path/data

module use /home/app/modulefiles/
module load python/cpu-3.7.4


python $loc/tools/get_everything_ready.py $sbid $path
bash $SCRIPTS/download_selavy.sh

ls $SCRIPTS/bash_GETDATA_beam*.sh | xargs -n 1 -P 8 --replace bash {}
ls $SCRIPTS/bash_GETDATA_beam*.sh | xargs -n 1 --replace bash {}

cd $DATA
ls *"$sbid"*tar | xargs -n 1 -P 10 --replace tar xvf {}

bash $loc/tools/submit_onebeam.sh beam$i $path 




