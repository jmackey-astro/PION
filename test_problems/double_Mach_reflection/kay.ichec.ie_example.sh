#!/bin/sh 

#SBATCH --time=00:20:00
# N.B. Kay has 20 processors per node, so 64-core job needs 4 nodes, etc.
#SBATCH --nodes=1 
#SBATCH -A dias01
#SBATCH -p DevQ
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jmackey@cp.dias.ie

module load intel

mkdir -p /ichec/work/dias01/jmackey/DMR
opdir=/ichec/work/dias01/jmackey/DMR

mpirun -np 16 ../../icgen-ug params_DMR_n520.txt silo
mpirun -np 16 ../../pion-ug DMRm10t60_n520_0000.00000000.silo \
  outfile=${opdir}/DMRm10t60_n520 redirect=${opdir}/log_DMRm10t60_n520

# NOTES:
# Submit job using:
#    sbatch kay.ichec.ie_example.sh
# Query the queue using:
#    squeue --account=dias01
# 

