
sbid=$1
beam=$2

path=/o9000/ASKAP/VAST/fast_survey/SB$sbid

bash $path/scripts/kill_beam"$beam"_jobs.sh
bash ~/vast_fastdetection/tools/clean_work_tree.sh $sbid $beam
bash $path/scripts/bash_GETDATA_beam"$beam".sh

cd $path/data
tar xvf scienceData.*.SB"$sbid".*.beam"$beam"_averaged_cal.leakage.ms.tar

cd ../../
bash ~/vast_fastdetection/tools/submit_onebeam.sh beam"$beam" $path
