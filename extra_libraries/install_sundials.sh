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

NCORES=4
export CC=gcc
export CXX=g++
export FC=gfortran
SHARED=YES
HDF5_LIBS="/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial"

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

COMPILE_SUNDIALS=yes


##################################
##########   SUNDIALS   ##########
##################################
if [ "$COMPILE_SUNDIALS" == "yes" ]
then
#################################
# Change these for new versions:
  VERSION=sundials-5.6.1
  FILE="${VERSION}".tar.gz
  SRC_DIR="${VERSION}"
  BLD_DIR=sundials_build
  REMOTE_URL=https://computing.llnl.gov/projects/sundials/download/"${FILE}"
  echo "********************************"
  echo "*** INSTALLING SUNDIALS/CVODE LIBRARY FILE=${FILE}****"
  echo "********************************"
#################################
  if [ -e $FILE ]; then
    echo "*** File downloaded, moving on... ***"
  else 
    echo "***** File does not exist ******"
    echo "********************************"
    echo "*** Automatic download of ${FILE}"
    echo "from ${REMOTE_URL}"
    $WGET --no-check-certificate $REMOTE_URL -O $FILE
    # check it downloaded.
    if [ ! -f $FILE ]; then
      echo "File not found! : $FILE"
      echo "Download of Silo Library Failed... quitting"
      exit
    fi
    echo "***  Downloaded."
    echo "********************************"
  fi 
  export CFLAGS='-O3'
  echo "***********************************"
  echo "*** EXTRACTING SUNDIALS LIBRARY ***"
  echo "***********************************"
  tar zxvf $FILE
  echo "***********************************"
  echo "*** RUNNING CMAKE CONFIG ***"
  echo "***********************************"
  BASE_PATH=`pwd`
  echo "Path = $BASE_PATH"
  mkdir -p $BLD_DIR
  cd $BLD_DIR
  echo "Running CMAKE"
  if [ "$SHARED" == "NO" ]
  then
    echo " ****** NOT COMPILING SHARED LIBRARIES ****** "
    cmake -DCMAKE_INSTALL_PREFIX=${BASE_PATH} \
   -DEXAMPLES_INSTALL_PATH=${BASE_PATH} -DEXAMPLES_INSTALL=OFF \
   -DBUILD_ARKODE=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF \
   -DBUILD_KINSOL=OFF -DBUILD_CVODES=OFF \
   -DBUILD_CVODE=ON -DBUILD_SHARED_LIBS=OFF \
   -DSUNDIALS_INSTALL_CMAKEDIR=share/sundials \
   ${BASE_PATH}/${SRC_DIR}
  else
    echo " ****** COMPILING SHARED LIBRARIES ****** "
    cmake -DCMAKE_INSTALL_PREFIX=${BASE_PATH} \
   -DEXAMPLES_INSTALL_PATH=${BASE_PATH} -DEXAMPLES_INSTALL=OFF \
   -DBUILD_ARKODE=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF \
   -DBUILD_KINSOL=OFF -DBUILD_CVODES=OFF \
   -DBUILD_CVODE=ON \
   -DSUNDIALS_INSTALL_CMAKEDIR=share/sundials \
   ${BASE_PATH}/${SRC_DIR}
  fi
  echo "********************************"
  echo "*** RUNNING MAKE ***"
  echo "********************************"
  make -j$NCORES
  make -j$NCORES install
  echo "*********************************"
  echo "*** INSTALLED CVODE LIBRARY ***"
  echo "*********************************"
  cd $CURDIR
  rm -rf $BLD_DIR
  echo "********************************"
  echo "*** FINISHED! ***"
  echo "********************************"
fi

cd $CURDIR

