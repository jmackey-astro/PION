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

COMPILE_SILO=yes
COMPILE_SUNDIALS=yes
COMPILE_FITS=yes

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
  COMPILE_SILO=yes
  COMPILE_SUNDIALS=yes
  COMPILE_FITS=no


elif [ "$id" == "Debian" ] && [ "$code" == "stretch" ]; then
  echo "Detected Debian 9 (stretch), use system libs for SILO, FITS, GSL, SUNDIALS!"
  echo "No need to compile extra libraries"
  echo "run  apt install libcfitsio-bin libcfitsio-dev libsilo-dev libsilo-bin python-silo libsundials-dev openmpi-bin openmpi-common curl libhdf5-serial-dev git libgsl-dev"
  echo "Then cd to serial/parallel binary directory and run  bash compile_code.sh"
  exit

elif [ "$id" == "Debian" ] && [ "$code" == "buster" ]; then
  echo "Detected Debian 10 (buster), use system libs for SILO, FITS, GSL, SUNDIALS!"
  echo "No need to compile extra libraries"
  echo "run  apt install libcfitsio-bin libcfitsio-dev libsilo-dev libsilo-bin python-silo libsundials-dev openmpi-bin openmpi-common curl libhdf5-serial-dev git libgsl-dev"
  echo "Then cd to serial/parallel binary directory and run  bash compile_code.sh"
  exit

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
  FILE=silo-4.10.2-bsd.tar.gz
  SRC_DIR=silo-4.10.2-bsd
  REMOTE_URL=https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz
#################################


  if [ -e $FILE ]; then
          echo "*** File exists, no need to download ***"
  else 
          echo "***** File does not exist ******"
          echo "********************************"
          echo "*** DOWNLOADING SILO LIBRARY ***"
          echo "********************************"
          if [ $MAKE_UNAME == "osx" ]; then
            curl  $REMOTE_URL -o $FILE
          else
            $WGET --no-check-certificate $REMOTE_URL -O $FILE
            if [ $? != 0 ]; then
              echo ** failed to download with wget, trying curl instead **
              rm $FILE
              curl $REMOTE_URL -o $FILE
            fi
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
#
  if [ "$MAKE_UNAME" == "KAY" ]
  then
    echo " ****** KAY.ICHEC.IE no shared libs, no python ****** "
    ./configure --prefix=${BASE_PATH} \
   --disable-browser \
   --disable-fortran \
   --disable-silex \
   --disable-shared \
   --disable-pythonmodule
  elif [ "$SHARED" == "NO" ]
  then
    echo " ****** NOT COMPILING SHARED LIBRARIES ****** "
    ./configure --prefix=${BASE_PATH} \
   --disable-browser \
   --disable-fortran \
   --disable-silex \
   --disable-shared \
   --enable-pythonmodule
  else
    echo " ****** COMPILING SHARED LIBRARIES ****** "
    ./configure --prefix=${BASE_PATH} \
   --disable-fortran \
   --disable-silex \
   --enable-pythonmodule \
   --with-hdf5=$HDF5_LIBS
  fi

  echo "********************************"
  echo "*** RUNNING MAKE ***"
  echo "********************************"
  make -j$NCORES
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
  bash install_sundials5.sh
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

