#!/usr/bin/env bash

beam=$1 # in format of beamxx

# Download and untar visibility
bash ~/vast_fastdetection/examples/bash_GETDATA_"$beam".sh

# Rescale data (from askapsoft convention to CASA convention)
# And fix beam positions
FIXDATA=`sbatch ~/vast_fastdetection/examples/slurm_FIXDATA_"$beam".sh | awk '{print $4}'`

# Create a deep sky model and subtract it
MODELING=`sbatch -d afterok:${FIXDATA} ~/vast_fastdetection/examples/slurm_MODELING_"$beam".sh | awk '{print $4}'`

# Generate short images 
IMGFAST=`sbatch -d afterok:${MODELING} ~/vast_fastdetection/examples/slurm_IMGFAST_"$beam".sh | awk '{print $4}'`

echo scancel $FIXDATA $MODELING $IMGFAST > /o9000/ASKAP/VAST/fast_survey/kill_"$beam"_jobs.sh
