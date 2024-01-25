#!/usr/bin/env bash

beam=$1 # in format of beamxx
path=$2

# move to SCRIPTS directory first 
# Download and untar visibility
#bash $path/scripts/bash_CHECKDATA.sh

# Rescale data (from askapsoft convention to CASA convention)
# And fix beam positions
FIXDATA=`sbatch "$path"/scripts/slurm_FIXDATA_"$beam".sh | awk '{print $4}'`

# Create a deep sky model and subtract it
#cd $path/models
#MODELING=`sbatch -d afterok:${FIXDATA} "$path"/scripts/slurm_MODELING_"$beam".sh | awk '{print $4}'`
#MODELING=`sbatch "$path"/scripts/slurm_MODELING_"$beam".sh | awk '{print $4}'`

# Generate short images 
cd $path/images
IMGFAST=`sbatch -d afterok:${FIXDATA} "$path"/scripts/slurm_IMGFAST_"$beam".sh | awk '{print $4}'`

# Candidates selection 
SELCAND=`sbatch -d afterok:${IMGFAST} "$path"/scripts/slurm_SELCAND_"$beam".sh | awk '{print $4}'`

# Clean intermediate data products
CLNDATA=`sbatch -d afterok:${SELCAND} "$path"/scripts/slurm_CLNDATA_"$beam".sh | awk '{print $4}'`


echo scancel $FIXDATA $IMGFAST $SELCAND $CLNDATA > $path/scripts/kill_"$beam"_jobs.sh
