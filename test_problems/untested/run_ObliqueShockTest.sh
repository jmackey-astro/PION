#!/bin/bash
#
# (c) Jonathan Mackey 2016.
#
# - 2016.05.21 JM: Run the Oblique Shock test, in parallel.

BASE_DIR=`pwd`
code_dir=${BASE_DIR}/../bin_parallel
test_dir=${BASE_DIR}
src_dir=${BASE_DIR}/../source
ICGEN=${BASE_DIR}/../icgen-ug
PION=${BASE_DIR}/../pion-ug



DATE=`date +%Y-%m-%d`

# aibn129 settings.
#data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/tests_$DATE
#visit_cmd=/vol/aibn129/aibn129_1/jmackey/extra_libraries/visit_bin/bin/visit
#visit_cmd=/vol/software/software/tools/visit/bin/visit

# OS-X settings.
#data_dir=/Users/jm/Documents/CODE/pion_dev/test_problems/data_$DATE
#visit_cmd=/Applications/VisIt.app/Contents/MacOS/VisIt

data_dir=${BASE_DIR}/data_${DATE}

mkdir -p $data_dir




#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Check we can compile the code
#
cd $code_dir
echo "MAKE IN" $code_dir
bash ./compile_code.sh
if [ ! -f ../pion-ug ] || [ ! -f ../icgen-ug ]
then
  echo "Cannot compile code"
  exit
else
  echo "MAKE SUCEEDED"
fi
#---------------------------------------------------------------------
#---------------------------------------------------------------------



#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Oblique Shock (Quirk Instability) Problem
#
echo "*************** Quirk Instability 2D Hydro ***************"
cd ${test_dir}/ObliqueShock/
./run_ObliqueShockTest.sh $test_dir $code_dir $data_dir ${ICGEN} ${PION}
#./make_OST_image.sh $test_dir $code_dir $data_dir $visit_cmd ObliqueM40_Roe_FKJav01 ObliqueM40_Roe_FKJav01
#./make_OST_image.sh $test_dir $code_dir $data_dir $visit_cmd ObliqueM40_Roe_Hcorr   ObliqueM40_Roe_Hcorr  
#./make_OST_image.sh $test_dir $code_dir $data_dir $visit_cmd ObliqueM25_Roe_FKJav01 ObliqueM25_Roe_FKJav01
#./make_OST_image.sh $test_dir $code_dir $data_dir $visit_cmd ObliqueM25_Roe_Hcorr   ObliqueM25_Roe_Hcorr  
#eog *.jpeg &
cd $test_dir
#
# Oblique Shock (Quirk Instability) Problem
#---------------------------------------------------------------------
#---------------------------------------------------------------------

exit












