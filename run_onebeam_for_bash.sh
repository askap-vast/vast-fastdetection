#!/bin/bash 

# bash run_onebeam_for_bash.sh 45965 02

sbid=$1
i=$2 # beam number 

path=`pwd`/SB$sbid
loc=$(dirname "$0")

echo "The script you are running has:"
echo "basename: [$(basename "$0")]"
echo "dirname : [$(dirname "$0")]"
echo "pwd     : [$(pwd)]"
echo
echo The output will store to $path
echo The scripts are located in $loc


SCRIPTS=$path/scripts

python $loc/tools/get_everything_ready_for_bash.py $sbid $path
bash $SCRIPTS/download_selavy.sh
bash $SCRIPTS/bash_PROCESSING_beam"$i".sh
echo $sbid $i Finished!! 



