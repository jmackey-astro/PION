#!/bin/bash
#
# Assumes:
# (1) That PION has been built in pion/build/ with SILO support, PARALLEL build
# (2) That the file "data/BW2d_StoneMHD_l2_B010_n128_level00_0000.00000000.silo"
#     exists for comparison
# (3) that "silocompare" has also been built, and is available to compare the two files

#if [[ $# -ne 1 ]]; then
#    echo -e "Please provide the Pion build directory"
#    exit 1
#fi

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"

mpirun --allow-run-as-root --oversubscribe -np 4 ../icgen-ng \
  ${script_dir}/problems/MHD_Blastwave2D/params_MHD_blastwave2D_l2_B010_n128.txt \
  silo redirect=iclog

mpirun --allow-run-as-root --oversubscribe -np 4 ../pion-ng \
  BW2d_StoneMHD_l2_B010_n128_level00_0000.00000000.silo solver=7 cfl=0.24 \
  opfreq_time=0.1 outfile=NG_B010_n128 redirect=pionlog

REF_FILE=NG_B010_n128_level00_0000.00000824.silo
NEW_FILE=`ls NG_B010_n128_level00_0000.*.silo | tail -n1`
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
