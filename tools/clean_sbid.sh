# clean all visibilities and images 

sbid=$1
num=1 # number of parallel processes 

path=`pwd`/SB$sbid
loc=$(dirname "$0")

echo "The script you are running has:"
echo "basename: [$(basename "$0")]"
echo "dirname : [$(dirname "$0")]"
echo "pwd     : [$(pwd)]"
echo The output will store to $path
echo The scripts are located in $loc

echo cleaning all visibilities
rm -r $path/data/scienceData*

echo cleaning all short images
rm $path/images/*fits

echo cleaning all model images
rm $path/models/*fits 

echo cleaning statistical images
rm $path/candidates/*{peak,std,chisquare,gaussian}.fits

echo Finished
echo 
