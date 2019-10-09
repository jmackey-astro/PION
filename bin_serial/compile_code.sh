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
#NCORES=1
export CXX=g++

#################################
## First try debian-based linux #
#################################
id=`lsb_release -s -i`
ver=`lsb_release -s -r`
code=`lsb_release -s -c`
nc=`nproc --all`

if [ "$id" == "Ubuntu" ] && [ "$ver" == "18.04" ]; then
  echo "Detected Ubuntu 18.04: Note system Silo library has a bug, you must install yourself using the script in PION/extra_libraries"
  MAKE_UNAME=ubuntu18
  export CXX=g++
  export CC=gcc
  export PION_OPTIONS="-DSERIAL -DSILO -DFITS"
  export PION_OPTIMISE=HIGH
  NCORES=$nc
elif [ "$id" == "Debian" ] && [ "$code" == "stretch" ]; then
  echo "Detected Debian 9 (stretch), using system libs for SILO, FITS, GSL, SUNDIALS"
  export PION_OPTIONS="-DSERIAL -DSILO -DFITS"
  export PION_OPTIMISE=HIGH
  #export PION_OPTIMISE=LOW
  #NCORES=1
  export CXX=g++
  NCORES=$nc
  MAKE_UNAME=debian9
elif [ "$id" == "Debian" ] && [ "$code" == "buster" ]; then
  echo "Detected Debian 10 (buster), using system libs for SILO, FITS, GSL, SUNDIALS"
  export PION_OPTIONS="-DSERIAL -DSILO -DFITS"
  export PION_OPTIMISE=HIGH
  #export PION_OPTIMISE=LOW
  #NCORES=1
  export CXX=g++
  NCORES=$nc
  MAKE_UNAME=debian10
else
  echo "Failed to find a known version of Linux: checking for other OS types."
  MAKE_UNAME=standard
fi
#################################

#################################
### TEST FOR OS X (DARWIN)    ###
#################################
DDD=`uname -a | grep "Darwin"`
if [ ! -z "$DDD" ]; then
  export PION_OPTIONS="-DSERIAL -DSILO -DFITS"
  export CXX=g++
  export CC=gcc
  echo "***** COMPILING WITH OS-X: host ${HOST}: COMPILERS ARE $CC $CXX "  
  MAKE_UNAME=OSX
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
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DINTEL"
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

#####################################################################
#### now compile the code:
#####################################################################
export MAKE_UNAME
echo "COMPILING WITH MACHINE: $MAKE_UNAME"
make -j${NCORES} uni
make -j${NCORES} NG
#####################################################################

exit 



