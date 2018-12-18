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

export WGET='wget'


##############################
## TEST FOR SuperMUC (2016) ##
##############################
case $HOST in
  login[0-9][0-9])
    echo "Compiling on SuperMUC"
    export CC=icc
    export CXX=icpc
    export FC=ifort
    MAKE_UNAME=SUPERMUC
    NCORES=8
    ;;
esac
#######################

##########################################
### TEST FOR DIRAC-2-COMPLEXITY (2015) ###
##########################################
case $HOSTNAME in
  dirac[0-9][0-9])
    echo "Compiling on DIRAC-Complexity"
    module list
    module load intel/compilers/13.0.0 intel/impi/4.1.0 intel/mkl/11.0.0 cfitsio
    MAKE_UNAME=DIRAC
    NCORES=8
    # -DINTEL means the code uses the intel math headers instead of gnu.
    export CC=icc
    export CXX=icpc
    export FC=ifort
  ;;
esac
#######################


##############################
### TEST FOR KAY.ICHEC.IE  ###
##############################
case $HOSTNAME in
  login[0-9].novalocal)
    echo "Compiling on KAY/ICHEC"
    source /usr/share/Modules/init/bash
    #module purge
    module load intel
    module load dev cmake3
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
 --enable-pythonmodule
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

##################################
##########   SUNDIALS   ##########
##################################

#################################
# Change these for new versions:
FILE=sundials-2.6.2.tar.gz
SRC_DIR=sundials-2.6.2
BLD_DIR=sundials_build
REMOTE_URL=https://computation.llnl.gov/projects/sundials/download/sundials-2.6.2.tar.gz
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

##################################
##########   CFITSIO    ##########
##################################
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
cd $CURDIR

