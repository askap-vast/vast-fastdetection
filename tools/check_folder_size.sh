#!/bin/sh

for i in SB*
do

  i=${i#SB}

  if [ -d SB"$i" ]
  then
    du -sh SB$i

  fi

done
