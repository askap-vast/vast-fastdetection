#!/bin/bash 

#SBATCH --partition=all-x86-cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=FIX01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodelist=hw-x86-cpu10
#SBATCH --mem=50gb
#SBATCH --output=/o9000/ASKAP/VAST/fast_survey/LOGS/slurm_FIXDATA_SB47253_beam01.output
#SBATCH --error=/o9000/ASKAP/VAST/fast_survey/LOGS/slurm_FIXDATA_SB47253_beam01.error
#SBATCH --export=all

module use /home/app/modulefiles
module load casacore/cpu-py3.6.5-3.1.0

time -p python ~/vast_fastdetection/tools/askapsoft_rescale.py /o9000/ASKAP/VAST/fast_survey/SB47253/visibilities/scienceData.VAST_1806-25.SB47253.VAST_1806-25.beam01_averaged_cal.leakage.ms /o9000/ASKAP/VAST/fast_survey/SB47253/visibilities/scienceData.VAST_1806-25.SB47253.VAST_1806-25.beam01_averaged_cal.leakage.ms.corrected

time -p python ~/vast_fastdetection/tools/fix_dir.py /o9000/ASKAP/VAST/fast_survey/SB47253/visibilities/scienceData.VAST_1806-25.SB47253.VAST_1806-25.beam01_averaged_cal.leakage.ms.corrected

