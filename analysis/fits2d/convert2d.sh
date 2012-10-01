#!/bin/bash
if [ "$1" = "" ] ||  [ "$2" = "" ] ||  [ "$3" = "" ] ||  [ "$4" = "" ]
then
 echo Usage: $0 filename step element-no file-description
 file=out_0001
 element=1
else
 file=$1
 step=$2
 element=$3
 desc=$4
fi
counter=0

make
for i in $( ls $file*.h5 ); 
do
src=$file.$counter.h5
dest=$file.$counter.$desc.fits
./h5tofits $src $dest $element
echo The counter is $counter
let counter=counter+$step
done
make clean
