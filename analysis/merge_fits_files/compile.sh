#!/bin/bash
if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ] || [ "$4" = "" ] || [ "$5" = "" ]
  then
  echo "Usage: $0 <infile-path> <infile-base> <outfilebase> <nproc> <make_uname>"
  echo "  <infile-path> can be a relative or absolute path to the directory containing the FITS files."
  echo "  <infile-base> is the FITS filename, minus the _0000.XXXXXXXX.fits part."
  echo "  <outfilebase> is the name of the merged fits files and should include a path, otherwise it will output to current directory."
  echo "  <nproc> is the number of MPI processes the simulation was run with, equal to the number of FITS files for each timestep."
  echo "  <make_uname> is the machine identifier: standard, OSX, other?"
  exit
else
  input_path=$1
  infile=$2
  outfile=$3
  nproc=$4
fi

MAKE_UNAME=$5

MAKE_UNAME=$MAKE_UNAME make -f Makefile -j4


#####################################################################
## fix some linking problem with OSX (this is new... 2016.05.25)
#####################################################################
echo "$MAKE_UNAME"
if [ "$MAKE_UNAME" = "OSX" ]; then
  echo "Fixing links"
  path=`pwd`
  install_name_tool -change libsundials_cvode.1.dylib      \
   ${path}/../../extra_libraries/lib/libsundials_cvode.1.dylib       \
   merge_fits_files
  install_name_tool -change libsundials_nvecserial.0.dylib \
   ${path}/../../extra_libraries/lib/libsundials_nvecserial.0.dylib \
   merge_fits_files 
fi
#####################################################################

exit


#valgrind  --leak-check=full --show-reachable=yes \
#./merge_fits_files $1 $2 $3 $4 $5 $6

#make clean

