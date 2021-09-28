#!/bin/bash
#
# 3D adiabatic hydro blastwave, with 3 levels of refinement.
# Checks that the hydro is working with NG, against a reference solution
# from master

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"

mpirun --allow-run-as-root --oversubscribe -np 4 ../icgen-ug \
  ${script_dir}/problems/blastwave_crt3d/params_BW_UGcrt3D_NR032.txt \
  silo redirect=iclog

mpirun --allow-run-as-root --oversubscribe -np 4 ../pion-ug \
  BW_UGcrt3D_NR032_0000.00000000.silo \
  outfile=BW_UGcrt3D_NR032_new redirect=pionlog


REF_FILE=BW_UGcrt3D_NR032_REF_0000.00000088.silo
NEW_FILE=`ls BW_UGcrt3D_NR032_new_0000.*.silo | tail -n1`
../silocompare . $NEW_FILE ${script_dir}/data $REF_FILE 0 cmp 2 > tmp.txt

if grep -q "RESULTS ARE THE SAME" tmp.txt; then
  echo -e "${GREEN}*** 3D-Blastwave TEST HAS BEEN PASSED ***"
  tail -n10 tmp.txt
  echo -e "*** 3D-Blastwave TEST HAS BEEN PASSED ***${NC}"
  exit 0
else
  echo -e "${RED}*** 3D-Blastwave TEST HAS BEEN FAILED ***"
  tail -n10 tmp.txt
  echo -e "*** 3D-Blastwave TEST HAS BEEN FAILED ***${NC}"
  exit 1 
fi
