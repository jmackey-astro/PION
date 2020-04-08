#!/bin/sh 

#SBATCH --time=00:20:00
# N.B. Kay has 40 processors per node, so 64-core job needs 2 nodes, etc.
#SBATCH --nodes=1 
#SBATCH -A dias01
#SBATCH -p DevQ
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jmackey@cp.dias.ie

module load intel

mkdir -p /ichec/work/dias01/jmackey/DTE2D
opdir=/ichec/work/dias01/jmackey/DTE2D

mpirun -np 16 ../../icgen-ug params_DTE_D2Full_TTI_n00256.txt silo
mpirun -np 32 ../../pion-ug DTE_D2Full_TTI_n00256_0000.00000000.silo \
  outfile=${opdir}/DTE_D2Full_TTI_n00256_s4 \
  redirect=${opdir}/log_DTE_D2Full_TTI_n00256_s4 \
  solver=4

# NOTES:
# Submit job using:
#    sbatch kay.ichec.ie_DTE2D_n00256.sh
# Query the queue using:
#    squeue --account=dias01
# 

