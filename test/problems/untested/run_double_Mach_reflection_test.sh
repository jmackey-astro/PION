#!/bin/bash
#
# Jonathan Mackey 2009.
# This script runs the double Mach reflection test for a given
# resolution.
#
# -2009-12-17 JM: Added in overstable shock test.
# -2009-12-23 JM: Added call to make shock tube plots.
# -2009-12-25 JM: Uncommented Shock tube tests.
# -2010-09-03 JM: updated directories for AIfA.
# -2010.10.11 JM: Updated directories to be relative to uniform_grid_code/
# -2010.12.06 JM: Added Field Loop test problem.
# -2012.07.06 JM: Modified to work in branches.
# -2016.05.24 SG: Modified to plot in Python instead of visit.


BASE_DIR=`pwd`
code_dir=${BASE_DIR}/../bin_serial
test_dir=${BASE_DIR}
src_dir=${BASE_DIR}/../source



#---------------------------------------------------------------------
# Default directories for storing data and finding the VisIt installation.
DATE=`date +%Y-%m-%d`
#visit_cmd='vglrun /usr/local/bin/visit'

# aibn129 settings.
#data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/tests_$DATE
#visit_cmd=/vol/software/software/tools/visit/bin/visit

# OS-X settings.
#data_dir=/Users/jm/Documents/CODE/pion_dev/test/problems/data_$DATE
#visit_cmd=/Applications/VisIt.app/Contents/MacOS/VisIt
#---------------------------------------------------------------------



#---------------------------------------------------------------------
# See if the required info can be parsed from the command-line.
# If not just guess.
#
if [ "$1" = "" ]
then
  echo "Please enter a resolution to use: NX = 130, 260, or 520"
  exit
else
  resolution=$1
fi

if [ "$2" = "" ]
then
  data_dir=${BASE_DIR}/data_$DATE
else
  data_dir=${BASE_DIR}/$2
fi

if [ "$3" = "" ]
then
  visit_cmd='vglrun /usr/local/bin/visit'
else
  visit_cmd=$3
fi
#---------------------------------------------------------------------

mkdir -p $data_dir
mkdir -p ${data_dir}/FIGS


#---------------------------------------------------------------------
# Set the code-testing flags correctly in testing_flags.h
#---------------------------------------------------------------------
#sed -i "" -e "s/^#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" ${src_dir}/defines/testing_flags.h
#sed -i "" -e "s/^#define CHECK_MAGP/\/\/#define CHECK_MAGP/g" ${src_dir}/defines/testing_flags.h
#sed -i "" -e "s/^#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/" ${src_dir}/defines/functionality_flags.h
#sed -i "" -e "s/^#define BLAST_WAVE_CHECK/\/\/#define BLAST_WAVE_CHECK/g" ${src_dir}/defines/testing_flags.h
#---------------------------------------------------------------------
# Run the Double Mach Reflection test
#---------------------------------------------------------------------
cd ${test_dir}/double_Mach_reflection
./run_DMR_tests.sh    ${test_dir}/double_Mach_reflection $code_dir ${data_dir}/DMR ${resolution}
./make_DMR_figures.sh ${test_dir}/double_Mach_reflection $code_dir ${data_dir}/DMR "$visit_cmd" ${resolution}

#---------------------------------------------------------------------
# Unset the TESTING flags so that the code is back to normal operation.
#sed -i -e "s/^#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" ${src_dir}/defines/testing_flags.h
#sed -i -e "s/^#define CHECK_MAGP/\/\/#define CHECK_MAGP/g" ${src_dir}/defines/testing_flags.h
#sed -i -e "s/^\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/" ${src_dir}/defines/functionality_flags.h
#sed -i -e "s/^#define BLAST_WAVE_CHECK/\/\/#define BLAST_WAVE_CHECK/g" ${src_dir}/defines/testing_flags.h
#---------------------------------------------------------------------

exit


