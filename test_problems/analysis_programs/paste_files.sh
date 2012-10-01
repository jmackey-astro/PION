#!/bin/bash
#

F1BASE=$1
F2BASE=$2

LIST1=(`ls ${F1BASE}*.txt`)
LIST2=(`ls ${F2BASE}*.txt`)

for ii in $(seq 0 $((${#LIST1[@]} - 1)))
do
  echo ${LIST1[$ii]} ${LIST2[$ii]} cmp_${LIST1[$ii]}
  paste ${LIST1[$ii]} ${LIST2[$ii]} > cmp_${LIST1[$ii]}
done

exit
