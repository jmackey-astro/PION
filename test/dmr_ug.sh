#!/bin/bash
#
# Author:  Jonathan Mackey
#
# 2021.05.14 JM: setting up a test calculation for pipeline
#
# Assumes:
# (1) That PION has been built in pion/build/ with SILO support, PARALLEL build
# (2) That the file "data/DMRm10t60_n160_0000.00000365.silo" exists for comparison
# (3) that "silocompare" has also been built, and is available to compare the two files

#if [[ $# -ne 1 ]]; then
#    echo "Please provide the Pion build directory"
#    exit 1
#fi

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"

mpirun --allow-run-as-root --oversubscribe -np 4 ../icgen-ug \
  ${script_dir}/problems/double_Mach_reflection/params_DMR_n160.txt silo
redirect=iclog |& grep -v "Read -1"

mpirun --allow-run-as-root --oversubscribe -np 4 ../pion-ug \
  DMRm10t60_n160_0000.00000000.silo outfile=DMR_new_n160 redirect=pionlog \
  |& grep -v "Read -1"


REF_FILE=DMRm10t60_n160_0000.00000365.silo
NEW_FILE=`ls DMR_new_n160_0000.*.silo | tail -n1`
../silocompare . $NEW_FILE ${script_dir}/data $REF_FILE 0 cmp 2 > tmp.txt

if grep -q "RESULTS ARE THE SAME" tmp.txt; then
  echo -e "${GREEN}*** TEST HAS BEEN PASSED ***"
  tail -n10 tmp.txt
  echo -e "*** TEST HAS BEEN PASSED ***${NC}"
  exit 0
else
  echo -e "${RED}*** TEST HAS BEEN FAILED ***"
  tail -n10 tmp.txt
  echo -e "*** TEST HAS BEEN FAILED ***${NC}"
  exit 1 
fi
