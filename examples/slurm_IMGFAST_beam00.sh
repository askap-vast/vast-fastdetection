#!/bin/bash 

#SBATCH --partition=all-x86-cpu
#SBATCH --time=30:00:00
#SBATCH --job-name=IMG-00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodelist=hw-x86-cpu10
#SBATCH --mem=30gb
#SBATCH --output=/o9000/ASKAP/VAST/fast_survey/LOGS/slurm_IMGFAST_SB47253_beam00.output
#SBATCH --error=/o9000/ASKAP/VAST/fast_survey/LOGS/slurm_IMGFAST_SB47253_beam00.error
#SBATCH --export=all

module use /home/app/modulefiles
module load casa/5.0.0-218.el6


cd /o9000/ASKAP/VAST/fast_survey/SB47253/images


time casa --logfile /o9000/ASKAP/VAST/fast_survey/LOGS/casa_IMGFAST_SB47253_beam00.log --nogui -c ~/vast_fastdetection/imaging/short_imaging.py /o9000/ASKAP/VAST/fast_survey/SB47253/visibilities/scienceData.VAST_1806-25.SB47253.VAST_1806-25.beam00_averaged_cal.leakage.ms.corrected SB47253_beam00 10

