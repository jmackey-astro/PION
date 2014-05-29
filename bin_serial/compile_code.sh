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
#
# We first need to set MAKE_UNAME which is an identifier for the computer
# we are compiling on.  If it is not a known computer we just set it to
# "standard" which is a standard linux system, assumed to have 4 cores (for
# compiling code only, not running).
#
MAKE_UNAME=standard
NCORES=4
export PION_OPTIONS="-DSERIAL -DSILO -DFITS"
export PION_OPTIMISE=HIGH
#export PION_OPTIMISE=LOW
export CXX=g++

#################################
### TEST FOR PHALANX ICC/ICPC ###
#################################
if [ "${HOST}" = 'phalanx.star.ucl.ac.uk' ]; then
  source /opt/intel/Compiler/11.1/073/bin/ifortvars.sh intel64
  source /opt/intel/Compiler/11.1/073/bin/iccvars.sh intel64
  export SGIMPT=/opt/sgi/mpt/mpt-1.26
  export PATH=$SGIMPT/bin:/opt/sgi/perfcatcher/bin:$PATH
  export LD_LIBRARY_PATH=$SGIMPT/lib:$LD_LIBRARY_PATH
  # -DINTEL means the code uses the intel math headers instead of gnu.
  export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DINTEL"
  export PION_OPTIMISE=HIGH
  export CXX=icpc
  echo "***** COMPILING WITH PHALANX: COMPILERS ARE $CC $CXX "  
  MAKE_UNAME=phalanx
  NCORES=8
fi
#################################

#################################
### TEST FOR Dougal ICC/ICPC ###
#################################
if [ "${HOST}" = 'dougal.hpc.phys.ucl.ac.uk' ]; then
  source /opt/intel/Compiler/11.1/046/bin/ifortvars.sh intel64
  source /opt/intel/Compiler/11.1/046/bin/iccvars.sh intel64
  # -DINTEL means the code uses the intel math headers instead of gnu.
  export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DINTEL"
  export PION_OPTIMISE=HIGH
  export CXX=icpc
  echo "***** COMPILING ON ${HOST}: COMPILERS ARE $CC $CXX "  
  MAKE_UNAME=dougal
  NCORES=8
fi
#################################

#######################
case $HOST in
  jj[0-9][0-9]l[0-9][0-9])
    echo "Compiling on JUROPA"
    module purge
    module load mkl/10.2.2.025 intel/11.1.059 sundials/2.4.0 parastation/mpi2-intel-5.0.25-2
    MAKE_UNAME=JUROPA
    NCORES=8
    # -DINTEL means the code uses the intel math headers instead of gnu.
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DINTEL"
    export PION_OPTIMISE=HIGH
    export CXX=icpc
    ;;
esac
#######################

#################################
### TEST FOR OS X (DARWIN)    ###
#################################
DDD=`uname -a | grep "Darwin"`
if [ ! -z "$DDD" ]; then
  export PION_OPTIONS="-DSERIAL -DSILO -DFITS"
  export PION_OPTIMISE=HIGH
  export CXX=g++
  export CC=gcc
  echo "***** COMPILING WITH OS-X: host ${HOST}: COMPILERS ARE $CC $CXX "  
  MAKE_UNAME=imac
  NCORES=2
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

#######################
### TEST FOR JUDGE ###
#######################
case $HOST in
  judgel[0-9])
    echo "Compiling on JUDGE"
    module purge
    module load intel/11.1.072 mkl/10.2.5.035 parastation/intel
    MAKE_UNAME=JUDGE
    NCORES=8
    # -DINTEL means the code uses the intel math headers instead of gnu.
    export PION_OPTIONS="-DSERIAL -DSILO -DFITS -DINTEL"
    export PION_OPTIMISE=HIGH
    export CXX=icpc
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

#####################################################################

#####################################################################
#### now compile the code:
#####################################################################
export MAKE_UNAME
echo "COMPILING WITH MACHINE: $MAKE_UNAME"
make -j${NCORES}
#####################################################################


exit



