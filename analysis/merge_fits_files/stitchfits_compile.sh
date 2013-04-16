#!/bin/bash
if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ] || [ "$4" = "" ]
  then
  echo "Usage: $0 <infile-path> <infile-base> <outfilebase> <nproc>"
  echo "  <infile-path> can be a relative or absolute path to the directory containing the FITS files."
  echo "  <infile-base> is the FITS filename, minus the _0000.XXXXXXXX.fits part."
  echo "  <outfilebase> is the name of the merged fits files and should include a path, otherwise it will output to current directory."
  echo "  <nproc> is the number of MPI processes the simulation was run with, equal to the number of FITS files for each timestep."
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

