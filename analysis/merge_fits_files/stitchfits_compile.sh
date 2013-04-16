#!/bin/bash
if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ] || [ "$4" = "" ]
  then
  echo Usage: $0 infile-path infile-base outfilebase nproc start-num step
  exit
else
  input_path=$1
  infile=$2
  outfile=$3
  nproc=$4
fi

make -f Makefile.stitchfits -j4

#valgrind  --leak-check=full --show-reachable=yes \
./stitchfits $1 $2 $3 $4 $5 $6

make -f Makefile.stitchfits clean

