#!/bin/bash
#
# Jonathan Mackey 2009.
# This script runs the 1D blastwave test.
#
# -2009-12-17 JM: Added in overstable shock test.
# -2009-12-23 JM: Added call to make shock tube plots.
# -2009-12-25 JM: Uncommented Shock tube tests.
# -2010-09-03 JM: updated directories for AIfA.
# -2010.10.11 JM: Updated directories to be relative to uniform_grid_code/
# -2010.12.06 JM: Added Field Loop test problem.
# -2012.07.06 JM: Modified to work in branches.
# - 2016.03.06 JM: updated for pion v1

BASE_DIR=`pwd`
code_dir=${BASE_DIR}/../bin_serial
test_dir=${BASE_DIR}
src_dir=${BASE_DIR}/../source


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
  data_dir=/Users/jm/Documents/CODE/pion_dev/test_problems/data_$DATE
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
#---------------------------------------------------------------------
# Test Problem: Blastwave (1D)
#
echo "*********** BLAST-WAVE TESTS!!! PUT YOUR GOGGLES ON ;-) **********"
cd ${test_dir}/blastwave_sph1d
sed -i -e "s|^\/\/\#define BLAST_WAVE_CHECK|\#define BLAST_WAVE_CHECK|g" \
  ${src_dir}/defines/testing_flags.h
##
## 1D spherical blast wave
##
./run_spherical_BW.sh $test_dir/blastwave_sph1d $code_dir ${data_dir}/blastwave_sph1d
echo "*********** 1D BLAST-WAVE TESTS FINISHED, MAKING FIGS NOW **********"
./make_spherical_BW_plots.sh ${data_dir}/blastwave_sph1d
echo "*********** 1D BW TESTS: FIGS GENERATED (SEE data-dir png/eps files) **********"
sed -i -e "s|^#define BLAST_WAVE_CHECK|//#define BLAST_WAVE_CHECK|g" ${src_dir}/defines/testing_flags.h
#
# End of Test Problem: Blastwave (1D)
#---------------------------------------------------------------------
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

exit




#exit







#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Oblique Shock (Quirk Instability) Problem
#
echo "*************** Quirk Instability 2D Hydro ***************"
cd ${test_dir}/ObliqueShock/
./run_ObliqueShockTest.sh $test_dir $code_dir $data_dir
./make_OST_image.sh $test_dir $code_dir $data_dir $visit_cmd ObliqueM40_Roe_FKJav01 ObliqueM40_Roe_FKJav01
./make_OST_image.sh $test_dir $code_dir $data_dir $visit_cmd ObliqueM40_Roe_Hcorr   ObliqueM40_Roe_Hcorr  
./make_OST_image.sh $test_dir $code_dir $data_dir $visit_cmd ObliqueM25_Roe_FKJav01 ObliqueM25_Roe_FKJav01
./make_OST_image.sh $test_dir $code_dir $data_dir $visit_cmd ObliqueM25_Roe_Hcorr   ObliqueM25_Roe_Hcorr  
eog *.jpeg &
cd $test_dir
#
# Oblique Shock (Quirk Instability) Problem
#---------------------------------------------------------------------
#---------------------------------------------------------------------



#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Test Problem:
#
#
# Overstable Shock 'Test'
#
cd ${test_dir}/test_OverstableShock
./run_RSH_tests.sh $test_dir $code_dir $data_dir
#
# Test Problem:
#---------------------------------------------------------------------
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# Tidy up stuff: move all figures to a directory in data_dir/
# Then make sure flags are unset in the code flags files.
#---------------------------------------------------------------------
cd ${test_dir}/
mkdir $data_dir/FIGS
mv *.jpeg *.eps */*.jpeg */*.eps $data_dir/FIGS

# Make sure these are unset at the end of the tests...
sed -i -e "s/^#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" \
${src_dir}/defines/testing_flags.h
sed -i -e "s/^#define CHECK_MAGP/\/\/#define CHECK_MAGP/g" \
${src_dir}/defines/testing_flags.h
sed -i -e "s/^#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_T/" \
${src_dir}/defines/functionality_flags.h
sed -i -e "s/^#define BLAST_WAVE_CHECK/\/\/#define BLAST_WAVE_CHECK/g" \
${src_dir}/defines/testing_flags.h

exit















#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Test Problem:
#

#
# Test Problem:
#---------------------------------------------------------------------
#---------------------------------------------------------------------


#####################################################################
########################## TEMPORARY STUFF ##########################
#####################################################################


# If the code is on a network drive, code may run faster with data_dir
# set to a local disk.
DATE=`date +%Y-%m-%d`
data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/tests_$DATE
mkdir $data_dir

# cmp_dir doesn't do anything yet.  It will contain some results from
# standard problems to compare the current tests to.
cmp_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests
#cmp_dir=${BASE_DIR}/trunk/test_results/ref

#visit_cmd=/vol/aibn129/aibn129_1/jmackey/extra_libraries/visit_bin/bin/visit
visit_cmd=/vol/software/software/tools/visit/bin/visit
#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Test Problem:
#
#
# Overstable Shock 'Test'
#
cd ${test_dir}/test_OverstableShock
./run_RSH_tests.sh $test_dir $code_dir $data_dir
#
# Test Problem:
#---------------------------------------------------------------------
#---------------------------------------------------------------------

exit
#####################################################################
########################## TEMPORARY STUFF ##########################
#####################################################################


