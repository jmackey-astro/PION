#!/bin/bash
#
# Jonathan Mackey 2016.
# This script runs the 3D blastwave test.
#
# - 2016.03.17 JM: wrote the script.
# - 2016.08.26 JM: updated to check consistency between serial and
#    parallel runs, and to plot the density vs. radius.

BASE_DIR=`pwd`
serial_dir=${BASE_DIR}/../bin_serial
pllel_dir=${BASE_DIR}/../bin_parallel
test_dir=${BASE_DIR}
src_dir=${BASE_DIR}/../source
cmp_dir=${BASE_DIR}/../analysis/silocompare


# Test command-line arguments
if test "$#" -ne 3; then
    echo "UASGE: run_blastwave_3d.sh <data_dir> <visit_command> <resolution [016/032/064/128/256]>"
    exit
fi

#---------------------------------------------------------------------
DATE=`date +%Y-%m-%d`
# OS-X settings.
data_dir=/Users/jm/Documents/CODE/pion_dev/test_problems/data_$DATE
visit_cmd=/Applications/VisIt.app/Contents/MacOS/VisIt
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# See if the directories are given from the command-line.

if [ "$1" = "" ]
then
  data_dir=./data_$DATE
else
  data_dir=$1
fi

if [ "$2" = "" ]
then
  visit_cmd=/Applications/VisIt.app/Contents/MacOS/VisIt
else
  visit_cmd=$2
fi

resolution=$3
#---------------------------------------------------------------------

mkdir -p $data_dir
mkdir -p ${data_dir}/FIGS

# use the correct version of sed
DDD=`uname -a | grep "Darwin"`
if [ ! -z "$DDD" ]; then
  sed_cmd=gsed
else
  sed_cmd=sed
fi


#---------------------------------------------------------------------
# Set the code-testing flags correctly in testing_flags.h
#---------------------------------------------------------------------
echo ${src_dir}
ls ${src_dir}/defines/testing_flags.h
${sed_cmd} -i -e "s/^#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" \
  ${src_dir}/defines/testing_flags.h
${sed_cmd} -i -e "s/^#define CHECK_MAGP/\/\/#define CHECK_MAGP/g" \
  ${src_dir}/defines/testing_flags.h
${sed_cmd} -i -e "s/^#define SET_NEGATIVE_PRES/\/\/#define SET_NEGATIVE_PRES/" \
  ${src_dir}/defines/functionality_flags.h
${sed_cmd} -i -e "s/^#define BLAST_WAVE_CHECK/\/\/#define BLAST_WAVE_CHECK/g" \
  ${src_dir}/defines/testing_flags.h

#---------------------------------------------------------------------
# Test Problem: Blastwave (3D)
#
echo "*********** BLAST--WAVE TESTS!!! PUT YOUR GOGGLES ON ;-) **********"
#
# Adiabatic 3D blast wave tests
#
# First run code in serial:
cd ${test_dir}/blastwave_crt3d
bash ./run_BW3D_serial.sh   $test_dir/blastwave_crt3d $serial_dir ${data_dir}/blastwave3D $resolution
bash ./run_BW3D_parallel.sh $test_dir/blastwave_crt3d $pllel_dir  ${data_dir}/blastwave3D $resolution
# compare the serial and parallel results.
bash ./compare_ser_pll.sh $test_dir/blastwave_crt3d ${cmp_dir}  ${data_dir}/blastwave3D $resolution
bash ./make_scatter_plots.sh ${data_dir}/blastwave3D $resolution

#
# End of Test Problem: Blastwave (3D)
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Unset the TESTING flags so that the code is back to normal operation.
#---------------------------------------------------------------------
${sed_cmd} -i -e "s/^#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" \
  ${src_dir}/defines/testing_flags.h
${sed_cmd} -i -e "s/^#define CHECK_MAGP/\/\/#define CHECK_MAGP/g" \
  ${src_dir}/defines/testing_flags.h
${sed_cmd} -i -e "s/^\/\/#define SET_NEGATIVE_PRES/#define SET_NEGATIVE_PRES/" \
  ${src_dir}/defines/functionality_flags.h
${sed_cmd} -i -e "s/^#define BLAST_WAVE_CHECK/\/\/#define BLAST_WAVE_CHECK/g" \
  ${src_dir}/defines/testing_flags.h
#---------------------------------------------------------------------






