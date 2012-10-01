#!/bin/bash
if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ] || [ "$4" = "" ] || [ "$5" = "" ]
then
echo Usage: $0 infilebase outfilebase nproc start-num step
else
infile=$1
outfile=$2
nproc=$3
start=$4
step=$5
fi
#cd ../testing/
#./clean.sh
cd ../analysis/
make -f Makefile.stitchfits
#valgrind  --leak-check=full --show-reachable=yes \
./stitchfits $1 $2 $3 $4 $5
make -f Makefile.stitchfits clean

