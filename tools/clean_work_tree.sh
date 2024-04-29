

# bash clean_work_tree.py 49156 20
sbid=$1
i=$2

#path=/o9000/ASKAP/VAST/fast_survey/SB$sbid
path=`pwd`/SB$sbid
echo cleaning work tree for SB$sbid beam$i in folder $path

rm -r $path/data/*beam$i*
rm -r $path/models/*beam$i*
rm -r $path/images/*beam$i*
rm -r $path/candidates/*beam$i*
