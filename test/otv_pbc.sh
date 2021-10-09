#!/bin/bash
#
# Author:  Jonathan Mackey
#
# 2021.08.13 JM: setting up a test calculation for pipeline.  This is
#  for the Orszag-Tang Vortex, with beta=10/5, Mach Number=1
#
# Assumes:
# (1) That PION has been built with SILO, PARALLEL, TOOLS, TEST-PROBLEMS support
# (2) That the file "" exists for comparison
# (3) that "silocompare" has also been built, and is available to compare the two files

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

###################################################################
# Logic for whether to test MPI, OpenMP, or hybrid parallelisation
# Default (also with no argument) is to use MPI with np=4
nt=0    # number of threads
mpi=""  # MPI command

if [[ $1 -eq 1 ]]; then
  nt=4
  mpi=""
elif [[ $1 -eq 2 ]]; then
  nt=2
  mpi="mpirun --allow-run-as-root --oversubscribe -np 2"
else
  nt=1
  mpi="mpirun --allow-run-as-root --oversubscribe -np 4"
fi
###################################################################

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"

${mpi} ../icgen-ug \
 ${script_dir}/problems/OrszagTangVortex/param_OrszagTang_n128.txt silo \
 redirect=iclog omp-nthreads=${nt}
${mpi} ../pion-ug \
 OrszagTang_n128_b3.33m1.0_0000.00000000.silo outfile=OrszagTang_n128_new \
 redirect=pionlog-otv opfreq_time=1.0 omp-nthreads=${nt}

REF_FILE=OrszagTang_n128_b3.33m1.0_0000.00001138.silo
NEW_FILE=`ls OrszagTang_n128_new_0000.*.silo | tail -n1`
../silocompare . $NEW_FILE ${script_dir}/data $REF_FILE 0 cmp 2 > tmp.txt

if grep -q "RESULTS ARE THE SAME" tmp.txt; then
  echo -e "${GREEN}*** OTV TEST HAS BEEN PASSED ***"
  tail -n10 tmp.txt
  echo -e "*** OTV TEST HAS BEEN PASSED ***${NC}"
  rm *.silo
  exit 0
else
  echo -e "${RED}*** OTV TEST HAS BEEN FAILED ***"
  tail -n10 tmp.txt
  echo -e "*** OTV TEST HAS BEEN FAILED ***${NC}"
  exit 1 
fi
