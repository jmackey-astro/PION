#!/bin/bash
#
# Jonathan Mackey 2009.
# This script runs the 2D blastwave test.
#
# -2009-12-17 JM: Added in overstable shock test.
# -2009-12-23 JM: Added call to make shock tube plots.
# -2009-12-25 JM: Uncommented Shock tube tests.
# -2010-09-03 JM: updated directories for AIfA.
# -2010.10.11 JM: Updated directories to be relative to uniform_grid_code/
# -2010.12.06 JM: Added Field Loop test problem.
# -2012.07.06 JM: Modified to work in branches.
# - 2016.03.14 JM: updated for pion v1

BASE_DIR=`pwd`
code_dir=${BASE_DIR}/../bin_serial
test_dir=${BASE_DIR}
src_dir=${BASE_DIR}/../source



#---------------------------------------------------------------------
DATE=`date +%Y-%m-%d`
# OS-X settings.
data_dir=/Users/jm/Documents/CODE/pion_dev/test/problems/data_$DATE
visit_cmd=/Applications/VisIt.app/Contents/MacOS/VisIt
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# See if the directories are given from the command-line.

if [ "$1" = "" ]
then
  data_dir=/Users/jm/Documents/CODE/pion_dev/test/problems/data_$DATE
else
  data_dir=$1
fi

if [ "$2" = "" ]
then
  visit_cmd=/Applications/VisIt.app/Contents/MacOS/VisIt
else
  visit_cmd=$2
fi
#---------------------------------------------------------------------

mkdir -p $data_dir
mkdir -p ${data_dir}/FIGS


#---------------------------------------------------------------------
# Set the code-testing flags correctly in testing_flags.h
#---------------------------------------------------------------------
sed -i "" -e "s/^#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" \
${src_dir}/defines/testing_flags.h
sed -i "" -e "s/^#define CHECK_MAGP/\/\/#define CHECK_MAGP/g" \
${src_dir}/defines/testing_flags.h
sed -i "" -e "s/^#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/" \
${src_dir}/defines/functionality_flags.h
sed -i "" -e "s/^#define BLAST_WAVE_CHECK/\/\/#define BLAST_WAVE_CHECK/g" \
${src_dir}/defines/testing_flags.h

#---------------------------------------------------------------------
# Test Problem: Blastwave (2D)
#
echo "*********** BLAST--WAVE TESTS!!! PUT YOUR GOGGLES ON ;-) **********"
#
# Adiabatic 2D blast wave tests
#
cd ${test_dir}/blastwave_axi2d
bash ./run_axisymmetric_BW.sh      $test_dir/blastwave_axi2d $code_dir ${data_dir}/blastwave_axi2d
#bash ./visit_BWaxi2d_make_plots.sh $test_dir/blastwave_axi2d $code_dir ${data_dir}/blastwave_axi2d $visit_cmd
#
# End of Test Problem: Blastwave (2D)
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Unset the TESTING flags so that the code is back to normal operation.
#---------------------------------------------------------------------
sed -i -e "s/^#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" \
${src_dir}/defines/testing_flags.h
sed -i -e "s/^#define CHECK_MAGP/\/\/#define CHECK_MAGP/g" \
${src_dir}/defines/testing_flags.h
sed -i -e "s/^\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/" \
${src_dir}/defines/functionality_flags.h
sed -i -e "s/^#define BLAST_WAVE_CHECK/\/\/#define BLAST_WAVE_CHECK/g" \
${src_dir}/defines/testing_flags.h
#---------------------------------------------------------------------






