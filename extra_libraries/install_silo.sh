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

export WGET='wget'

case $HOSTNAME in
  login[0-9].kay.ichec.ie)
    echo "Compiling on KAY/ICHEC"
    #source /usr/share/Modules/init/bash
    #module purge
######### gcc #########
#    module load cmake3
#    module load gcc
######### gcc #########
######## intel ########
    #module load cmake3
    #module load intel
    #module load gcc/8.2.0
    #export CC=icc
    #export CXX=icpc
    #export FC=ifort
######## intel ########
    HDF5=NO
    SHARED=YES
    module list
    ;;
esac

DDD=`uname -a | grep "Darwin"`
if [ ! -z "$DDD" ]; then
  #export CXX=g++-12
  #export CC=gcc-12
  echo "***** COMPILING WITH OS-X: host ${HOST}: COMPILERS ARE $CC $CXX "  
  MAKE_UNAME=OSX
  #NCORES=1
  path=`pwd`
  export PYTHON=/usr/local/bin/python3
fi


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

  ./configure --prefix=${BASE_PATH} \
    --disable-browser \
    --disable-fortran \
    --disable-silex \
    --disable-pythonmodule --enable-shared

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

