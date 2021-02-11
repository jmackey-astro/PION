#!/bin/bash
#
# - 2012.02.22 JM/HD: Added options for compiling on UCL-dougal
# - 2012.09.11 JM: Added options for SuperMUC
# - 2013.02.07 JM: Updated for new library versions (all).
# - 2014.04.14 JM: Added section for JUDGE at JSC.
# - 2015.01.14 JM: Section for Juropatest system at JSC.
# - 2016.04.29 JM: updated for sundials 2.6.2 and cfitsio 3390.
# - 2016.05.04 JM: Added FIONN to list of machines.
# - 2016.05.25 JM: Added support for OSX (use curl to download).

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
#module load cmake3/3.12.3
#module load python py/intel
#module load python numpy
module list

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
  FILE=silo-4.10.2-bsd.tgz
  SRC_DIR=silo-4.10.2-bsd
  REMOTE_URL=https://wci.llnl.gov/sites/wci/files/2021-01/silo-4.10.2-bsd.tgz
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
    --disable-pythonmodule

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
