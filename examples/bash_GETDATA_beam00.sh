#!/bin/bash

# Generate automatically from a python script
# Download and untar visibility for SB47253 beam00
# You can run this in terminal directly, simply use "bash examples/bash_GETDATA_beam00.sh" 



echo 
wget -O scienceData.VAST_1806-25.SB47253.VAST_1806-25.beam00_averaged_cal.leakage.ms.tar https://casda.csiro.au/casda_data_access/data/sync?id=DTrJQv6PDYh7jGflwu_qcntEhjTzCJChqZHtZFc53Y_u65sO82TKoseLxV7HiG9P -t 0
wget -O scienceData.VAST_1806-25.SB47253.VAST_1806-25.beam00_averaged_cal.leakage.ms.tar https://casda.csiro.au/casda_data_access/data/sync?id=DTrJQv6PDYh7jGflwu_qcntEhjTzCJChqZHtZFc53Y_u65sO82TKoseLxV7HiG9P -t 0 -c

echo 
tar xvf scienceData.VAST_1806-25.SB47253.VAST_1806-25.beam00_averaged_cal.leakage.ms.tar

