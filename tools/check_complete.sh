#!/bin/sh

for i in SB*
do

  i=${i#SB}

  if [ -d SB"$i" ]
  then
    num=$(ls SB"$i"/candidates/*peak_cand.csv 2> /dev/null | wc -l)
    cand=$(ls SB"$i"/candidates/*lightcurve*png 2> /dev/null | wc -l)
    echo SB$i beam=$num cands=$cand

    if (( $num != 36 ))
    then
      echo WARNING: SB$i $num is not completed!!!
    fi

  fi

done
