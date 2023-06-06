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

ml OpenMPI/4.1.1-GCC-10.3.0
ml git
ml HDF5/1.10.7-gompi-2021a
ml spdlog/1.9.2-GCCcore-10.3.0
ml CMake/3.20.1-GCCcore-10.3.0
ml Boost/1.76.0-GCC-10.3.0

NCORES=4
export CC=mpicc
export CXX=mpicxx
export FC=mpifort
SHARED=YES
HDF5=YES
HDF5_LIBS="/apps/all/HDF5/1.10.7-gompi-2021a/include,/apps/all/HDF5/1.10.7-gompi-2021a/lib"

COMPILE_SILO=yes

export WGET='wget'
export NCORES
CURDIR=`pwd`

##################################
##########     SILO     ##########
##################################
if [ "$COMPILE_SILO" == "yes" ]
then

  echo "********************************"
  echo "*** INSTALLING SILO LIBRARY ****"
  echo "********************************"
#################################
# Change these for new versions:
  FILE=silo-4.11-bsd-smalltest.tar.gz
  SRC_DIR=silo-4.11-bsd
  REMOTE_URL=https://github.com/LLNL/Silo/releases/download/v4.11/silo-4.11-bsd-smalltest.tar.gz
#################################


  if [ -e $FILE ]; then
          echo "*** File exists, no need to download ***"
  else 
          echo "***** File does not exist ******"
          echo "********************************"
          echo "*** DOWNLOADING SILO LIBRARY ***"
          echo "********************************"
          $WGET --no-check-certificate $REMOTE_URL -O $FILE
          if [ $? != 0 ]; then
            echo ** failed to download with wget, trying curl instead **
            rm $FILE
            curl $REMOTE_URL -o $FILE
          fi
          # check it downloaded.
          if [ ! -f $FILE ]; then
            echo "File not found! : $FILE"
            echo "Download of Silo Library Failed... quitting"
            exit
          fi
  fi 
#echo "********************************"
#echo "*** EXTRACTING SILO LIBRARY ***"
#echo "********************************"
  tar zxf $FILE
#echo "********************************"
#echo "*** RUNNING CONFIGURE ***"
#echo "********************************"
  BASE_PATH=`pwd`
  echo "***Path = $BASE_PATH ***"
  cd $SRC_DIR
  make clean

  ./configure --prefix=/home/it4i-jmackey/code/libs \
    --disable-browser \
    --disable-fortran \
    --disable-silex \
    --disable-pythonmodule --enable-shared \
    --with-hdf5=$HDF5_LIBS

  echo "********************************"
  echo "*** RUNNING MAKE ***"
  echo "********************************"
  make -j $NCORES
  echo "********************************"
  echo "*** INSTALLING SILO LIBRARY ***"
  echo "********************************"
  make install
  cd $CURDIR
  echo "********************************"
  echo "*** FINISHED! ***"
  echo "********************************"
fi


cd $CURDIR

