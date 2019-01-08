#!/bin/sh 

#SBATCH --time=01:00:00
# N.B. Kay has 40 processors per node, so 64-core job needs 2 nodes, etc.
#SBATCH --nodes=1 
#SBATCH -A dias01
#SBATCH -p DevQ
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jmackey@cp.dias.ie

module load intel

mkdir -p /ichec/work/dias01/jmackey/DTE2D
opdir=/ichec/work/dias01/jmackey/DTE2D

mpirun -np 32 ../../icgen_NG_parallel params_DTE_NG_D2_TTI_n00128.txt silo
mpirun -np 32 ../../pion_NG_parallel DTE_NG_D2_TTI_n00128_level00_0000.00000000.silo \
  outfile=${opdir}/DTE_NG_D2_TTI_n00128 opfreq=512 \
  redirect=${opdir}/log_DTE_NG_D2_TTI_n00128_s6 solver=6

# NOTES:
# Submit job using:
#    sbatch kay.ichec.ie_DTE2D_n00128.sh
# Query the queue using:
#    squeue --account=dias01
# 

