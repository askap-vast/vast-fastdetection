
sbid=$1

path=/o9000/ASKAP/VAST/fast_survey/SB$sbid

for i in {00..35}
  do
    bash $path/scripts/kill_beam"$i"_jobs.sh
  done


rm -r $path

echo Finished
