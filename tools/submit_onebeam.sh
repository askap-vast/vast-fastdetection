#!/usr/bin/env bash

beam=$1 # in format of beamxx

# move to SCRIPTS directory first 
# Download and untar visibility
bash bash_GETDATA_"$beam".sh

# Rescale data (from askapsoft convention to CASA convention)
# And fix beam positions
FIXDATA=`sbatch slurm_FIXDATA_"$beam".sh | awk '{print $4}'`

# Create a deep sky model and subtract it
MODELING=`sbatch -d afterok:${FIXDATA} slurm_MODELING_"$beam".sh | awk '{print $4}'`

# Generate short images 
IMGFAST=`sbatch -d afterok:${MODELING} slurm_IMGFAST_"$beam".sh | awk '{print $4}'`

echo scancel $FIXDATA $MODELING $IMGFAST > kill_"$beam"_jobs.sh
