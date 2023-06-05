#!/bin/bash


#SBATCH --partition=all-x86-cpu
#SBATCH --time=10:00:00
#SBATCH --job-name=SEL-00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodelist=purley-x86-cpu05
#SBATCH --mem=50gb
#SBATCH --output=/o9000/ASKAP/VAST/fast_survey/LOGS/slurm_SELCAND.output
#SBATCH --error=/o9000/ASKAP/VAST/fast_survey/LOGS/slurm_SELCAND.error
#SBATCH --export=all

# setting python profile
module use /home/app/modulefiles
module load python/cpu-3.7.4

time python ~/vast_fastdetection/run_all.py /o9000/ASKAP/VAST/fast_survey/SB47253/models/SB47253_beam33.image.tt0.fits /o9000/ASKAP/VAST/fast_survey/SB47253/data/selavy-image.i.VAST_1806-25.SB47253.cont.taylor.0.restored.conv.components.xml /o9000/ASKAP/VAST/fast_survey/SB47253/images/ beam33 /o9000/ASKAP/VAST/fast_survey/SB47253/candidates/ SB47253_beam33 





