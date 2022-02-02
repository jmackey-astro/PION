#!/bin/bash

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

####################################
# Logic for whether to test MPI, OpenMP, or hybrid parallelisation
nt=0    # number of threads
mpi=""  # MPI command

if [[ $1 -eq 1 ]]; then
  nt=4
  mpi=""
elif [[ $1 -eq 2 ]]; then
  nt=2
  mpi="mpirun --allow-run-as-root --oversubscribe -np 2"
elif [[ -z $1 ]]; then
  nt=1
  mpi="mpirun --allow-run-as-root --oversubscribe -np 4"
else
  nt=1
  mpi="mpirun --allow-run-as-root --oversubscribe -np $1"
fi

echo " **** DTE-TEST: arg=$1  nt=$nt  mpi=$mpi **** "
####################################

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"

${mpi} ../icgen-ug  ${script_dir}/problems/DTE2D/params_DTEHD_d2l1n0064.txt \
    silo omp-nthreads=${nt} || exit 1
${mpi} ../pion-ug DTEHD_d2l1n0064_0000.00000000.silo outfile=dtehd_d2l1n0064 \
    finishtime=1.5e13 opfreq_time=0.5e13 omp-nthreads=${nt} || exit 1

REF_FILE=dtehd_d2l1n0064_0000.00002653.silo
NEW_FILE=`ls dtehd_d2l1n0064_0000.*.silo | tail -n1`
../silocompare . $NEW_FILE ${script_dir}/data $REF_FILE 0 cmp 2 > tmp.txt

if grep -q "RESULTS ARE THE SAME" tmp.txt; then
  echo "${GREEN}*** D-type HII region TEST HAS BEEN PASSED ***"
  tail -n10 tmp.txt
  echo "*** D-type HII region TEST HAS BEEN PASSED ***${NC}"
  exit 0
else
  echo "${RED}*** D-type HII region TEST HAS BEEN FAILED ***"
  tail -n10 tmp.txt
  echo "*** D-type HII region TEST HAS BEEN FAILED ***${NC}"
  exit 1 
fi
