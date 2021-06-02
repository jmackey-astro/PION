#!/bin/bash
#
# Author:  Jonathan Mackey
#
# 2021.05.14 JM: setting up a test calculation for pipeline
#
# Assumes:
# (1) That PION has been built in pion/build/ with SILO support, PARALLEL build
# (2) That the file "data/DMRm10t60_n160_0000.00000365.silo" exists for comparison
# (2) that "silocompare" has also been built, and is available to compare the two files

#if [[ $# -ne 1 ]]; then
#    echo "Please provide the Pion build directory"
#    exit 1
#fi

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"

../icgen-ug ${script_dir}/problems/double_Mach_reflection/params_DMR_n160.txt silo redirect=iclog
../pion-ug DMRm10t60_n160_0000.00000000.silo outfile=DMR_new_n160 redirect=pionlog

REF_FILE=DMRm10t60_n160_0000.00000365.silo
NEW_FILE=`ls DMR_new_n160_0000.*.silo | tail -n1`
../silocompare . $NEW_FILE ${script_dir}/data $REF_FILE 0 cmp 2 > tmp.txt

if grep -q "RESULTS ARE THE SAME" tmp.txt; then
  echo "*** TEST HAS BEEN PASSED ***"
  tail -n10 tmp.txt
  rm tmp.txt *.silo iclog_0_info.txt pionlog_0_info.txt cmp_0.txt
  echo "*** TEST HAS BEEN PASSED ***"
  exit 0
else
  echo "*** TEST HAS BEEN FAILED ***"
  tail -n10 tmp.txt
  echo "*** TEST HAS BEEN FAILED ***"
  #rm tmp.txt *.silo
  exit 1 
fi
