#!/bin/bash
#
# 3D adiabatic hydro blastwave, with uniform grid.
# The blast is off-centre, and we have periodic boundaries, so the shockwave
# crosses the periodic boundaries.
# Checks that the hydro is working with UG, against a reference solution
# from master

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
  ${script_dir}/problems/blastwave_crt3d/params_BW_UGcrt3Dpbc_NR032.txt \
  silo omp-nthreads=${nt} || exit 1

${mpi} ../pion-ug \
  BW_UGcrt3Dpbc_NR032_0000.00000000.silo \
  outfile=BW_UGcrt3Dpbc_NR032_new omp-nthreads=${nt} || exit 1


REF_FILE=BW_UGcrt3Dpbc_NR032_REF_0000.00000140.silo
NEW_FILE=`ls BW_UGcrt3Dpbc_NR032_new_0000.*.silo | tail -n1`
../silocompare . $NEW_FILE ${script_dir}/data $REF_FILE 0 cmp 2 > tmp.txt

if grep -q "RESULTS ARE THE SAME" tmp.txt; then
  echo -e "${GREEN}*** 3Dpbc-Blastwave TEST HAS BEEN PASSED ***"
  tail -n10 tmp.txt
  echo -e "*** 3Dpbc-Blastwave TEST HAS BEEN PASSED ***${NC}"
  exit 0
else
  echo -e "${RED}*** 3Dpbc-Blastwave TEST HAS BEEN FAILED ***"
  tail -n10 tmp.txt
  echo -e "*** 3Dpbc-Blastwave TEST HAS BEEN FAILED ***${NC}"
  exit 1 
fi
