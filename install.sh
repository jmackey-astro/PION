#!/bin/bash

# INSTALL PION
# Author: Jonathan Mackey
# Date: 2011-2016

BASE_DIR=`pwd`

echo "************** INSTALLING LIBRARIES **************"
cd ${BASE_DIR}/extra_libraries/
rm -rf bin lib include
bash ./install_all_libs.sh
echo "************** ----------------------------- *************"
echo "************** FINISHED INSTALLING LIBRARIES **************"
echo "************** ----------------------------- *************"
echo "COMPILING PARALLEL CODE"
cd ${BASE_DIR}/bin_parallel
bash ./clean.sh
bash ./compile_code.sh

#echo "COMPILING SERIAL CODE"
#cd ${BASE_DIR}/bin_serial
#bash ./clean.sh
#bash ./compile_code.sh
#echo "FINISHED COMPILING CODE"

cd ${BASE_DIR}
echo "Serial and parallel executables should be in current directory:"
echo " $BASE_DIR "
echo "If, not, then run the steps in install.sh one-by-one to see \
what the problem is."

exit



