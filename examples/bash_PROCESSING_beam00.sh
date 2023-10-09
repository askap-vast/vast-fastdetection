#!/bin/bash

# Generate automatically from a python script
# Prcessing data for SB49585 beam00
# It will download calibrated visibilities from CASDA
# You can run this in terminal directly, simply use "bash /import/eve/fast_imaging/SB49585/scripts/bash_PROCESSING_beam00.sh" 



echo beam00: Download calibrated visibilities data
wget -O /import/eve/fast_imaging/SB49585/data/scienceData.VAST_1510-56.SB49585.VAST_1510-56.beam00_averaged_cal.leakage.ms.tar https://casda.csiro.au/casda_data_access/data/sync?id=33gpOl7Qbz6Ir4TtEZeOXNt_bMveA7C1d5YPCTjPIPvaANahlWwqPS7JcVcsGgjP -t 0 -c

echo beam00: Untar the data to measurement sets
tar xvf /import/eve/fast_imaging/SB49585/data/scienceData.VAST_1510-56.SB49585.VAST_1510-56.beam00_averaged_cal.leakage.ms.tar -C /import/eve/fast_imaging/SB49585/data

echo beam00: Fix the measurement sets flux scaling
python /home/ywan3191/scripts/vast_fastdetection/tools/askapsoft_rescale.py /import/eve/fast_imaging/SB49585/data/scienceData.VAST_1510-56.SB49585.VAST_1510-56.beam00_averaged_cal.leakage.ms /import/eve/fast_imaging/SB49585/data/scienceData.VAST_1510-56.SB49585.VAST_1510-56.beam00_averaged_cal.leakage.ms.corrected

echo beam00: Fix the measurement sets pointing
python /home/ywan3191/scripts/vast_fastdetection/tools/fix_dir.py /import/eve/fast_imaging/SB49585/data/scienceData.VAST_1510-56.SB49585.VAST_1510-56.beam00_averaged_cal.leakage.ms.corrected

echo beam00: Create sky model and subtract...
cd /import/eve/fast_imaging/SB49585/models 
casa --log2term --logfile /import/eve/fast_imaging/SB49585/logfiles/casa_MODELING_SB49585_beam00.log --nogui -c /home/ywan3191/scripts/vast_fastdetection/imaging/model_making.py /import/eve/fast_imaging/SB49585/data/scienceData.VAST_1510-56.SB49585.VAST_1510-56.beam00_averaged_cal.leakage.ms.corrected SB49585_beam00

echo beam00: Create model-subtracted short images...
cd /import/eve/fast_imaging/SB49585/images 
casa --log2term --logfile /import/eve/fast_imaging/SB49585/logfiles/casa_IMGFAST_SB49585_beam00.log --nogui -c /home/ywan3191/scripts/vast_fastdetection/imaging/short_imaging.py /import/eve/fast_imaging/SB49585/data/scienceData.VAST_1510-56.SB49585.VAST_1510-56.beam00_averaged_cal.leakage.ms.corrected SB49585_beam00 10

echo beam00: Select candidates...
python /home/ywan3191/scripts/vast_fastdetection/run_all.py /import/eve/fast_imaging/SB49585/models/SB49585_beam00.image.tt0.fits /import/eve/fast_imaging/SB49585/data/selavy-image.i.VAST_1510-56.SB49585.cont.taylor.0.restored.conv.components.xml /import/eve/fast_imaging/SB49585/images beam00 /import/eve/fast_imaging/SB49585/candidates SB49585_beam00

echo beam00: Clean intermediate products...
# rm /import/eve/fast_imaging/SB49585/data/scienceData.VAST_1510-56.SB49585.VAST_1510-56.beam00_averaged_cal.leakage.ms.tar
rm -r /import/eve/fast_imaging/SB49585/data/scienceData.VAST_1510-56.SB49585.VAST_1510-56.beam00_averaged_cal.leakage.ms
mv /import/eve/fast_imaging/SB49585/models/*beam00*.fits /import/eve/fast_imaging/SB49585/fitsfiles
mv /import/eve/fast_imaging/SB49585/images/*beam00*.fits /import/eve/fast_imaging/SB49585/fitsfiles
rm -r /import/eve/fast_imaging/SB49585/models/*beam00*
rm -r /import/eve/fast_imaging/SB49585/images/*beam00*
mv /import/eve/fast_imaging/SB49585/fitsfiles/*beam00*image*fits /import/eve/fast_imaging/SB49585/models
mv /import/eve/fast_imaging/SB49585/fitsfiles/*beam00*fits /import/eve/fast_imaging/SB49585/images

echo beam00: Finished!