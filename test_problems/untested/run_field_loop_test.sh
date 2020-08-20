#!/bin/bash
#
# Jonathan Mackey 2009.
# This script runs all of my standard test problems in a series of
# subdirectories.  More tests will be added as I have time to work
# on it.
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
# aibn129 settings.
data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/tests_$DATE
visit_cmd=/vol/software/software/tools/visit/bin/visit
# OS-X settings.
data_dir=/Users/jm/Documents/CODE/pion_dev/test_problems/data_$DATE
visit_cmd=/Applications/VisIt.app/Contents/MacOS/VisIt
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# See if the directories are given from the command-line.
if [ "$1" = "" ]
then
  echo "Please enter a resolution to use: NX = 100, 200"
  exit
else
  resolution=$1
fi

if [ "$2" = "" ]
then
  data_dir=/Users/jm/Documents/CODE/pion_dev/test_problems/data_$DATE
else
  data_dir=$2
fi

if [ "$3" = "" ]
then
  visit_cmd=/Applications/VisIt.app/Contents/MacOS/VisIt
else
  visit_cmd=$3
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
# Run the Field Loop test
#---------------------------------------------------------------------
echo "*************** FIELD LOOP TEST 2D MHD ***************"
cd ${test_dir}/FieldLoop/
#
# We set a flag so that the code outputs the magnetic pressure on the
# full domain as a function of time.
#
sed -i -e "s|^\/\/#define CHECK_MAGP|#define CHECK_MAGP|g" \
  ${src_dir}/defines/testing_flags.h
#
# Now run the tests:
#
./run_FL_test.sh ${test_dir}/FieldLoop $code_dir $data_dir $visit_cmd ${resolution}
#
# Unset the special FieldLoop flag
#
sed -i -e "s|^#define CHECK_MAGP|\/\/#define CHECK_MAGP|g" \
  ${src_dir}/defines/testing_flags.h
echo "*************** FIELD LOOP TEST DONE   ***************"
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




#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Run the Ray-tracing tests (no dynamics)
#
# We want to set a testing flag, so the code does not calculate
# boundary-point timesteps in the first 10 timesteps, and also sets a
# fixed cross-section for hydrogen ionisation and recombination.
#
sed -i -e "s|//#define RT_TEST_PROBS|#define RT_TEST_PROBS|g" ${src_dir}/defines/testing_flags.h
#cd ${test_dir}/RT_1D_norec/
#./run_A1A3_tests.sh    $test_dir $code_dir $data_dir $visit_cmd
#cd ${test_dir}/RT_1D_rec/
#./run_A1A3_tests.sh    $test_dir $code_dir $data_dir $visit_cmd
#cd ${test_dir}/test_RT_nodyn/
#./run_RTnodyn_tests.sh $test_dir $code_dir $data_dir $visit_cmd
# Unset testing flag.
sed -i -e "s|^#define RT_TEST_PROBS|//#define RT_TEST_PROBS|g" ${src_dir}/defines/testing_flags.h
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Run the Double Mach Reflection test
#---------------------------------------------------------------------
cd ${test_dir}/test_DMR/
./run_DMR_tests.sh    $test_dir $code_dir $data_dir
./make_DMR_figures.sh $test_dir $code_dir $data_dir $visit_cmd
#---------------------------------------------------------------------
#---------------------------------------------------------------------





#---------------------------------------------------------------------
#---------------------------------------------------------------------
# 1D/2D shock tubes
#
echo "************** RUNNING  1D/2D SHOCK-TUBES **************"
cd ${test_dir}/test_ShockTubes/
# WARNING: HI-RES SHOCK-TUBE TESTS TAKE HOURS TO FINISH.
bash ./run_hiresST_tests.sh    $test_dir $code_dir $data_dir/STtest
bash ./run_ST2D_tests.sh $test_dir $code_dir $data_dir/STtest
bash ./run_ST_tests.sh   $test_dir $code_dir $data_dir/STtest
bash ./ST2D_analyse_tests.sh $test_dir $code_dir $data_dir/STtest
bash ./make_ST_plots.sh  ${test_dir}/test_ShockTubes $code_dir $data_dir/STtest
echo "************** FINISHED 1D/2D SHOCK-TUBES **************"
#---------------------------------------------------------------------
#---------------------------------------------------------------------



#exit

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Test Problem: Blastwave (1D,2D (slab+axi),3D)
#
echo "*********** BLAST--WAVE TESTS!!! PUT YOUR GOGGLES ON ;-) **********"
cd ${test_dir}/blastwave_sph1d
sed -i -e "s|//#define BLAST_WAVE_CHECK|#define BLAST_WAVE_CHECK|g" ${src_dir}/defines/testing_flags.h
##
## 1D spherical blast wave
##
./run_spherical_BW.sh $test_dir/blastwave_sph1d $code_dir ${data_dir}/blastwave_sph1d
echo "*********** 1D BLAST--WAVE TESTS FINISHED, MAKING FIGS NOW **********"
./make_spherical_BW_plots.sh ${data_dir}/blastwave_sph1d
mv $test_dir/blastwave_sph1d/*.jpeg $test_dir
echo "*********** 1D BW TESTS: FIGS GENERATED (SEE blastwave_sph1d/*.jpeg) **********"
sed -i -e "s|^#define BLAST_WAVE_CHECK|//#define BLAST_WAVE_CHECK|g" ${src_dir}/defines/testing_flags.h
#
# Adiabatic 2D blast wave tests
#
cd ${test_dir}/blastwave_axi2d
./run_axisymmetric_BW.sh      $test_dir/blastwave_axi2d $code_dir ${data_dir}/blastwave_axi2d
./visit_BWaxi2d_make_plots.sh $test_dir/blastwave_axi2d $code_dir ${data_dir}/blastwave_axi2d $visit_cmd
#
# End of Test Problem: Blastwave (1D,2D (slab+axi),3D)
#---------------------------------------------------------------------
#---------------------------------------------------------------------







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


