#!/bin/bash
#
# - 2021.04.20 JM: script to install latest boost libs.

mkdir include
mkdir bin
mkdir lib

NCORES=4
export CC=gcc
export CXX=g++
export FC=gfortran
SHARED=YES

export WGET='wget'

case $HOSTNAME in
  login[0-9].kay.ichec.ie)
    echo "Compiling on KAY/ICHEC"
    source /usr/share/Modules/init/bash
    module purge
    module load cmake3
    module load gcc
#module load cmake3/3.12.3
#module load python py/intel
#module load python numpy
    module list
    ;;
esac


export NCORES
CURDIR=`pwd`

COMPILE_BOOST=yes


##################################
##########   BOOST   ##########
##################################
if [ "$COMPILE_BOOST" == "yes" ]
then
#################################
# Change these for new versions:
  FILE=boost_1_75_0.tar.bz2
  VERSION=1.75.0
  BLD_DIR=boost_build
  REMOTE_URL=https://dl.bintray.com/boostorg/release/1.75.0/source/${FILE}
  echo "********************************"
  echo "*** INSTALLING BOOST FILE=${FILE}****"
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
      echo "Download of Boost Library Failed... quitting"
      exit
    fi
    echo "***  Downloaded."
    echo "********************************"
  fi 
  export CFLAGS='-O3'
  echo "***********************************"
  echo "*** EXTRACTING BOOST LIBRARY ***"
  echo "***********************************"
  tar --bzip2 -xf $FILE
  BASE_PATH=`pwd`
  echo "Path = $BASE_PATH"
  cd boost_1_75_0
  #./bootstrap.sh --show-libraries
  ./bootstrap.sh --prefix=${BASE_PATH}/boost --with-libraries=math
  ./b2 install
  echo "********************************"
  echo "*** FINISHED! ***"
  echo "********************************"
fi

cd $CURDIR

