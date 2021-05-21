#!/bin/bash
#
# Authors: Jonathan Mackey, Harpreet Dhanoa, others?
#
# - 2012-2020: ongoing development to add options for new machines.
# 

mkdir include
mkdir bin
mkdir lib
mkdir lib64
mkdir share
mkdir cmake

NCORES=4
export CC=gcc
export CXX=g++
export FC=gfortran
SHARED=YES
HDF5=YES
HDF5_LIBS="/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial"

COMPILE_SILO=yes
COMPILE_SUNDIALS=yes
COMPILE_FITS=yes

export WGET='wget'

source /usr/share/Modules/init/bash
#module purge
module load cmake3
module load gcc
    module load conda
    source activate
    module list
    MAKE_UNAME=KAY
    NCORES=8
    export CC=icc
    export CXX=icpc
    export FC=ifort
    SHARED=NO
    . ./install_python_silo.sh
    COMPILE_SILO=yes
    COMPILE_SUNDIALS=yes
    COMPILE_FITS=no

export NCORES
CURDIR=`pwd`


##################################
##########     BOOST    ##########
##################################
bash install_boost.sh

##################################
##########     SILO     ##########
##################################
if [ "$COMPILE_SILO" == "yes" ]
then
  bash install_silo.sh
fi


##################################
##########   SUNDIALS   ##########
##################################
if [ "$COMPILE_SUNDIALS" == "yes" ]
then
  bash install_sundials.sh
fi

##################################
##########   CFITSIO    ##########
##################################
if [ "$COMPILE_FITS" == "yes" ]
then
  echo "*******************************"
  echo "*** INSTALLING FITS LIBRARY ***"
  echo "*******************************"
#################################
# Change these for new versions:
  FILE=cfitsio3390.tar.gz
  SRC_DIR=cfitsio
  REMOTE_URL=https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio3390.tar.gz
#################################

  if [ -e $FILE ]; then
          echo "*** File exists, no need to download ***"
  else 
    echo "***** File does not exist ******"
    echo "*******************************"
    echo "*** DOWNLOADING FITS LIBRARY ***"
    echo "*******************************"
    $WGET --no-check-certificate $REMOTE_URL -O $FILE
    # check it downloaded.
    if [ ! -f $FILE ]; then
      echo "File not found! : $FILE"
      echo "Download of Silo Library Failed... quitting"
      exit
    fi
  fi
  echo "*******************************"
  echo "*** EXTRACTING FITS LIBRARY ***"
  echo "*******************************"
  tar zxf $FILE
  echo "*******************************"
  echo "*** RUNNING CONFIGURE ***"
  echo "*******************************"
  BASE_PATH=`pwd`
  echo "***Path = $BASE_PATH ***"
  cd $SRC_DIR
  make clean
  ./configure --prefix=${BASE_PATH}
  echo "*******************************"
  echo "*** RUNNING MAKE ***"
  echo "*******************************"
  make -j$NCORES
  echo "*******************************"
  echo "*** INSTALLING CFITSIO ***"
  echo "*******************************"
  make install
  echo "*******************************"
  echo "*** FINISHED! ***"
  echo "*******************************"
fi

cd $CURDIR
