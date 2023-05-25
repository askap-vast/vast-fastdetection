#!/bin/bash 

#SBATCH --partition=all-x86-cpu
#SBATCH --time=30:00:00
#SBATCH --job-name=MOD01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodelist=hw-x86-cpu10
#SBATCH --mem=200gb
#SBATCH --output=/o9000/ASKAP/VAST/fast_survey/LOGS/slurm_MODELING_SB47253_beam01.output
#SBATCH --error=/o9000/ASKAP/VAST/fast_survey/LOGS/slurm_MODELING_SB47253_beam01.error
#SBATCH --export=all

module use /home/app/modulefiles
module load casa/5.0.0-218.el6
module load python/cpu-3.6.5

# mkdir /o9000/ASKAP/VAST/fast_survey/SB47253/models
cd /o9000/ASKAP/VAST/fast_survey/SB47253/models

time casa --logfile /o9000/ASKAP/VAST/fast_survey/LOGS/casa_MODELING_SB47253_beam01.log --nogui -c ~/vast_fastdetection/imaging/model_making.py /o9000/ASKAP/VAST/fast_survey/SB47253/visibilities/scienceData.VAST_1806-25.SB47253.VAST_1806-25.beam01_averaged_cal.leakage.ms.corrected SB47253_beam01

