#!/bin/bash
#
# 2012.08.06 JM: Modified so it runs all sims at the same time.
#
cd ../../bin_serial
./compile_code.sh
cp ../bin/main_serial ../bin/icgen_serial ../test_problems/Raymond79_shock/
cd ../test_problems/Raymond79_shock/
#rm RSH1D_n*.*
#
./icgen_serial pf_RShRay79_wB_HHe_ModelE_n128.txt silo
./icgen_serial pf_RShRay79_wB_HHe_ModelE_n256.txt silo
./icgen_serial pf_RShRay79_wB_HHe_ModelE_n512.txt silo
./icgen_serial pf_RShRay79_wB_HHe_ModelE_n1024.txt silo
#
./main_serial IC_RSH1D_n128_v100_Ray79E_HHe.silo 5 1 optype=text \
 outfile=./RSH1D_n128_v100_Ray79E_HHe_C4 cooling=4 op_criterion=1 \
 opfreq_time=3.16e9 finishtime=3.16e11 > msg_RSH1D_n128_v100_Ray79E_HHe.txt &
#
./main_serial IC_RSH1D_n256_v100_Ray79E_HHe.silo 5 1 optype=text \
 outfile=./RSH1D_n256_v100_Ray79E_HHe_C4 cooling=4 op_criterion=1 \
 opfreq_time=3.16e9 finishtime=3.16e11 > msg_RSH1D_n256_v100_Ray79E_HHe.txt &
#
./main_serial IC_RSH1D_n512_v100_Ray79E_HHe.silo 5 1 optype=text \
 outfile=./RSH1D_n512_v100_Ray79E_HHe_C4 cooling=4 op_criterion=1 \
 opfreq_time=3.16e9 finishtime=3.16e11 > msg_RSH1D_n512_v100_Ray79E_HHe.txt &
#
./main_serial IC_RSH1D_n1024_v100_Ray79E_HHe.silo 5 1 optype=text \
 outfile=./RSH1D_n1024_v100_Ray79E_HHe_C4 cooling=4 op_criterion=1 \
 opfreq_time=1.58e9 finishtime=3.16e11 > msg_RSH1D_n1024_v100_Ray79E_HHe.txt &
#
wait
exit

