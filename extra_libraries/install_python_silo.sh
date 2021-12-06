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


COMPILE_SILO=yes
COMPILE_SUNDIALS=no
COMPILE_FITS=no

export WGET='wget'

#################################
## First try debian-based linux #
#################################
id=`lsb_release -s -i`
ver=`lsb_release -s -r`
code=`lsb_release -s -c`
nc=`nproc --all`

#################################

##############################
### TEST FOR KAY.ICHEC.IE  ###
##############################
case $HOSTNAME in
  login[0-9].kay.ichec.ie)
    echo "Compiling on KAY/ICHEC"
    #source /usr/share/Modules/init/bash
######### gcc #########
#    module load cmake3
#    module load gcc
######### gcc #########
######## intel ########
    #module load cmake3
    #module load intel
    #export CC=icc
    #export CXX=icpc
    #export FC=ifort
######## intel ########
    module load conda
    source activate
    module list
    MAKE_UNAME=KAY
    NCORES=8
    SHARED=YES
    COMPILE_SILO=yes
    COMPILE_SUNDIALS=no
    COMPILE_FITS=no
    ;;
esac
#######################

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
  FILE=silo-4.11-bsd.tgz
  SRC_DIR=silo-4.11-bsd
  REMOTE_URL=https://wci.llnl.gov/sites/wci/files/2021-09/silo-4.11-bsd.tgz
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
  BASE_PATH=${BASE_PATH}/python
  mkdir -p $BASE_PATH
  echo "***Path = $BASE_PATH ***"
  cd $SRC_DIR
  make clean
#
  if [ "$MAKE_UNAME" == "KAY" ]
  then
    echo " ****** KAY.ICHEC.IE shared python lib ****** "
    ./configure --prefix=${BASE_PATH} \
   --enable-browser \
   --disable-fortran \
   --disable-silex \
   --enable-shared \
   --enable-pythonmodule
  fi

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

