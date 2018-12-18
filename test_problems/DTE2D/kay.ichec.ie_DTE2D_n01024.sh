#!/bin/sh 

#SBATCH --time=04:00:00
# N.B. Kay has 20 processors per node, so 64-core job needs 4 nodes, etc.
#SBATCH --nodes=7
#SBATCH -A dias01
#SBATCH -p ProdQ
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jmackey@cp.dias.ie

module load intel

mkdir -p /ichec/work/dias01/jmackey/DTE2D
opdir=/ichec/work/dias01/jmackey/DTE2D

mpirun -np 16 ../../icgen_parallel params_DTE_D2Full_TTI_n01024.txt silo
mpirun -np 128 ../../pion_parallel DTE_D2Full_TTI_n01024_0000.00000000.silo \
  outfile=${opdir}/DTE_D2Full_TTI_n01024_s4 \
  redirect=${opdir}/log_DTE_D2Full_TTI_n01024_s4 \
  solver=4

# NOTES:
# Submit job using:
#    sbatch kay.ichec.ie_DTE2D_n01024.sh
# Query the queue using:
#    squeue --account=dias01
# 

