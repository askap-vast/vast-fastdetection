#!/usr/bin/env bash

# Download and untar visibility
bash ~/vast_fastdetection/examples/bash_GETDATA_beam00.sh

# Rescale data (from askapsoft convention to CASA convention)
# And fix beam positions
FIXDATA=`sbatch ~/vast_fastdetection/examples/slurm_FIXDATA.sh | awk '{print $4}'`

# Create a deep sky model and subtract it
MODELING=`sbatch -d afterok:${FIXDATA} ~/vast_fastdetection/examples/slurm_MODELING.sh | awk '{print $4}'`

echo scancel $FIXDATA $MODELING > /o9000/ASKAP/VAST/fast_survey/kill_beam00_jobs.sh
