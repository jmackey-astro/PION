#!/bin/bash
#
# Author: Jonathan Mackey
# Date:   2010.XX.XX
# Description:
#  Compilation script for the serial version of the code.
# 
# Modifications:
# - 2010.11.09 JM: Added options for Phalanx at UCL
# - 2011.12.20 JM: Changed completely to work with new Makefile.
#     This script now makes the decisions about compilers and optimisation.
# - 2012.02.22 JM/HD: Added options for Dougal at UCL
# - 2012.03.29 JM: Added options for OS-X.
# - 2012.09.11 JM: Added options for SuperMUC
# - 2013.01.17 JM: Got rid of readline/ncurses from link line in
#    production version of pion.
# - 2013.02.15 JM: Added NEW_METALLICITY flag for testing the new
#    microphysics classes.
# - 2013.02.27 JM: Added extensions for contributed code.
# - 2014.04.14 JM: Added option for JUDGE.
# - 2016.05.04 JM: Added FIONN to list of machines
# - 2020.11.25 JM: Updated regularly for new machines and settings.
#
# We first need to set MAKE_UNAME which is an identifier for the computer
# we are compiling on.  If it is not a known computer we just set it to
# "standard" which is a standard linux system, assumed to have 8 cores (for
# compiling code only, not running).
#
MAKE_UNAME=standard
#MAKE_UNAME=locallibs
export PION_OPTIONS="-DSERIAL -DSILO -DFITS"
export PION_OPTIMISE=HIGH
NCORES=8
#export PION_OPTIMISE=LOW
NCORES=1
export CXX=g++

#################################
## First try debian-based linux #
#################################
LINUX=""
id=""
ver=""
code=""
nc=""
if ! [ -x "$(command -v lsb_release)" ]; then
  LINUX="NO"
else
  LINUX="YES"
  id=`lsb_release -s -i`
  ver=`lsb_release -s -r`
  code=`lsb_release -s -c`
  nc=`nproc --all`
fi

if [ "$LINUX" == "YES" ]; then

  if [ "$id" == "Ubuntu" ] && [ "$ver" == "20.04" ]; then
    echo "Detected Ubuntu 20.04: Note system Silo library has a bug, you must install yourself using the script in PION/extra_libraries"
    MAKE_UNAME=debian10
    export CXX=g++
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DCVODE3"
    export PION_OPTIMISE=HIGH
    NCORES=$nc
    #NCORES=1
  elif [ "$id" == "Ubuntu" ] && [ "$ver" == "18.04" ]; then
    echo "Detected Ubuntu 18.04 (bionic): Note system Silo library has a bug, you must install yourself using the script in PION/extra_libraries"
    MAKE_UNAME=ubuntu18
    export CXX=g++
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DCVODE5"
    export PION_OPTIMISE=HIGH
    #export PION_OPTIMISE=LOW
    NCORES=$nc
  elif [  "$id" == "Ubuntu" ] && [ "$ver" == "16.04" ]; then
    echo "Detected Ubuntu 16.04 (xenial): compiling extra libraries"
    MAKE_UNAME=ubuntu16
    export CXX=g++
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DCVODE5"
    export PION_OPTIMISE=HIGH
    NCORES=$nc
  elif [ "$id" == "Debian" ] && [ "$code" == "stretch" ]; then
    echo "Detected Debian 9 (stretch), using system libs for SILO, FITS, GSL, SUNDIALS"
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DCVODE2"
    export PION_OPTIMISE=HIGH
    #export PION_OPTIMISE=LOW
    #NCORES=1
    export CXX=g++
    NCORES=$nc
    MAKE_UNAME=debian9
  elif [ "$id" == "Debian" ] && [ "$code" == "buster" ]; then
    echo "Detected Debian 10 (buster), using system libs for SILO, FITS, GSL, SUNDIALS"
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DCVODE3"
    export PION_OPTIMISE=HIGH
    #export PION_OPTIMISE=LOW
    #NCORES=1
    NCORES=$nc
    export CXX=g++
    MAKE_UNAME=debian10
  else
    echo "Failed to find a known version of Linux: checking for other OS types."
    MAKE_UNAME=standard
  fi

fi
#################################

#################################
### TEST FOR OS X (DARWIN)    ###
#################################
DDD=`uname -a | grep "Darwin"`
if [ ! -z "$DDD" ]; then
  if ! [ -x "$(command -v /opt/local/bin/port)" ]; then
    SL="Homebrew"
  else
    SL="MacPorts"
  fi
  echo "detected support software ${SL} and assuming that is what you are using"
  if [ "$SL" == "MacPorts" ]; then
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DCVODE3"
    export CXX=g++
    export CC=gcc
    export PION_OPTIMISE=HIGH
    echo "*** COMPILING WITH OS-X: host ${HOST}: libs=${SL}: COMPILERS ARE $CC $CXX"
    echo "Make sure you installed silo, gsl, sundials, cfitsio with MacPorts"
    MAKE_UNAME=OSX-MP
  elif [ "$SL" == "Homebrew" ]; then
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DCVODE5"
    export CXX=g++
    export CC=gcc
    export PION_OPTIMISE=HIGH
    echo "*** COMPILING WITH OS-X: host ${HOST}: libs=${SL}: COMPILERS ARE $CC $CXX"
    echo "Make sure you installed gsl, sundials, cfitsio with Homebrew, and SILO in ../extra_libraries/"
    MAKE_UNAME=OSX-HB
  else
    echo "ERROR: Need Macports or Homebrew installed on OSX"
    exit
  fi
  NCORES=`sysctl -n hw.ncpu`
fi
#################################

#########################
### TEST FOR SuperMUC ###
#########################
case $HOST in
  login[0-9][0-9])
    echo "Compiling on SuperMUC"
    echo "MKL LIBS: $MKL_SHLIB"
    export CC=icc
    export CXX=icpc
    export FC=if90
    MAKE_UNAME=SUPERMUC
    NCORES=8
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DINTEL -DCVODE2"
    export PION_OPTIMISE=HIGH
  ;;
esac
#######################

#####################################################################
# For testing/debugging, we need to add -DTESTING to the compile   ##
# flags, and also to add readline and maybe ncurses to the linker. ##
#####################################################################
if [ $PION_OPTIMISE == LOW ]
  then
  echo "LDFLAGS= $DLFLAGS"
  export LDFLAGS=" -lreadline -lncurses "
  echo "LDFLAGS= $LDFLAGS"
  export PION_OPTIONS="$PION_OPTIONS -DTESTING"
else
  export LDFLAGS=""
fi
#####################################################################

#####################################################################
############ EXTRA CODE EXTENSIONS FROM CONTRIBUTED CODE ############
#####################################################################
# Harpreet Dhanoa's chemistry/microphysics module
#PION_OPTIONS+=" -DHARPREETS_CODE_EXT"

# Read in turbulence simulations provided by Blakesley Burkhart
#PION_OPTIONS+=" -DBBTURBULENCE_CODE_EXT"
#export PION_OPTIONS
#echo PION_OPTIONS: $PION_OPTIONS

#PION_OPTIONS+=" -DCODE_EXT_HHE"
#export PION_OPTIONS
#echo PION_OPTIONS: $PION_OPTIONS

PION_OPTIONS+=" -DCODE_EXT_SBII"
export PION_OPTIONS
echo PION_OPTIONS: $PION_OPTIONS

#PION_OPTIONS="$PION_OPTIONS -DLEGACY_CODE"
#export PION_OPTIONS
#####################################################################

# if we failed to detect the system, try assuming self-compiled libs
# are present.
if [ $MAKE_UNAME == "" ]; then
  MAKE_UNAME=locallibs
fi

#####################################################################
#### now compile the code:
#####################################################################
export MAKE_UNAME
echo "COMPILING WITH MACHINE: $MAKE_UNAME"
make -j${NCORES} uni
make -j${NCORES} NG
#####################################################################

exit 



