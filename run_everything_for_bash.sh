#!/bin/bash 

sbid=$1
num=18 # number of parallel processes 

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
#bash $SCRIPTS/bash_CHECKDATA.sh

ls $SCRIPTS/bash_PROCESSING_beam*.sh | xargs -n 1 -t -P $num bash

echo $sbid Finished!! 

# uplode final file and email alert
python $loc/get_overall_table.py $sbid
python $loc/tools/email_alert.py $sbid -e yuanmingwang@swin.edu.au

# run dyspec for known sources
python $loc/tools/run_dyspec_sbid.py $sbid

