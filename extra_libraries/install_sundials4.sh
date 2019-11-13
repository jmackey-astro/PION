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

MAKE_UNAME=standard
NCORES=4
export CC=gcc
export CXX=g++
export FC=gfortran
SHARED=YES
HDF5_LIBS="/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial"

export WGET='wget'

#################################
## First try debian-based linux #
#################################
id=`lsb_release -s -i`
ver=`lsb_release -s -r`
code=`lsb_release -s -c`
nc=`nproc --all`

if [ "$id" == "Ubuntu" ] && [ "$ver" == "18.04" ]; then
  echo "Detected Ubuntu 18.04: compiling extra libraries"
  MAKE_UNAME=ubuntu18
  export CXX=g++
  export CC=gcc
  export PION_OPTIONS="-DSERIAL -DSILO -DFITS"
  export PION_OPTIMISE=HIGH
  NCORES=$nc
  export CC=gcc
  export CXX=g++
  export FC=gfortran
  SHARED=YES
  HDF5_LIBS="/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial"

elif [ "$id" == "Debian" ] && [ "$code" == "stretch" ]; then
  echo "Detected Debian 9 (stretch)"

elif [ "$id" == "Debian" ] && [ "$code" == "buster" ]; then
  echo "Detected Debian 10 (buster)"

else
  echo "Failed to find a known version of Linux: checking for other OS types."
fi
#################################

##############################
### TEST FOR KAY.ICHEC.IE  ###
##############################
case $HOSTNAME in
  login[0-9].kay.ichec.ie)
    echo "Compiling on KAY/ICHEC"
    source /usr/share/Modules/init/bash
    #module purge
    module load cmake3
    module load intel/2018u4
    #module load cmake3/3.12.3
    #module load python py/intel
    #module load python numpy
    module list
    MAKE_UNAME=KAY
    NCORES=8
    export CC=icc
    export CXX=icpc
    export FC=ifort
    SHARED=NO
    ;;
esac
#######################

#################################
### TEST FOR OS X (DARWIN)    ###
#################################
DDD=`uname -a | grep "Darwin"`
if [ ! -z "$DDD" ]; then
  export CXX=g++
  export CC=gcc
  export FC=gfortran
  echo "***** COMPILING WITH OS-X: host ${HOST}: COMPILERS ARE $CC $CXX "  
  MAKE_UNAME=osx
  NCORES=4
  SHARED=NO
fi
#################################


export MAKE_UNAME
export NCORES
echo "COMPILING WITH MACHINE: ${MAKE_UNAME}. Compilers: CC=$CC FC=$FC CXX=$CXX"
CURDIR=`pwd`

COMPILE_SUNDIALS=yes


##################################
##########   SUNDIALS   ##########
##################################
if [ "$COMPILE_SUNDIALS" == "yes" ]
then
#################################
# Change these for new versions:
  FILE=sundials-4.1.0.tar.gz
  SRC_DIR=sundials-4.1.0
  BLD_DIR=sundials_build
  REMOTE_URL=https://computing.llnl.gov/projects/sundials/download/sundials-4.1.0.tar.gz
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
    if [ $MAKE_UNAME == "osx" ]; then
      curl  $REMOTE_URL -o $FILE
    else
      $WGET --no-check-certificate $REMOTE_URL -O $FILE
    fi
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
   ${BASE_PATH}/${SRC_DIR}
  else
    echo " ****** COMPILING SHARED LIBRARIES ****** "
    cmake -DCMAKE_INSTALL_PREFIX=${BASE_PATH} \
   -DEXAMPLES_INSTALL_PATH=${BASE_PATH} -DEXAMPLES_INSTALL=OFF \
   -DBUILD_ARKODE=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF \
   -DBUILD_KINSOL=OFF -DBUILD_CVODES=OFF \
   -DBUILD_CVODE=ON \
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

