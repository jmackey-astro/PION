#!/bin/bash

make -f Makefile.serial.code; make -f Makefile.serial.icgenerator
##################################
#### Stellar winds 2010-07-24 ####
##################################
../bin/icgen_serial pf_test_winds.txt silo
../bin/icgen_serial pf_test_winds2.txt silo
../bin/icgen_serial pf_test_winds3.txt silo
../bin/icgen_serial pf_test_winds4.txt silo

#valgrind --log-file=wind1_test_log 
../bin/pion_serial IC_wind_test1.silo 5 1 op_criterion=1 opfreq=1.58e10 cooling=0 redirect=test1 outfile=/export/aibn214_1/jmackey/testing/results/wind_test1 cfl=0.1 &
#valgrind --log-file=wind2_test_log 
../bin/pion_serial IC_wind_test2.silo 5 1 op_criterion=1 opfreq=1.58e10 cooling=0 redirect=test2 outfile=/export/aibn214_1/jmackey/testing/results/wind_test2 cfl=0.1 &
wait
#valgrind --log-file=wind3_test_log 
../bin/pion_serial IC_wind_test3.silo 5 1 op_criterion=1 opfreq=1.58e10 cooling=0 redirect=test3 outfile=/export/aibn214_1/jmackey/testing/results/wind_test3 cfl=0.1 &
#valgrind --log-file=wind4_test_log 
../bin/pion_serial IC_wind_test4.silo 5 1 op_criterion=1 opfreq=1.58e10 cooling=0 redirect=test4 outfile=/export/aibn214_1/jmackey/testing/results/wind_test4 cfl=0.1 &
wait
exit



##################################
## Overstable Radiative Shocks  ##
## 2009-12-16                   ##
##################################
../bin/icgen_serial test_RSH/pf_rsh_v140.txt silo
../bin/icgen_serial test_RSH/pf_rsh_chem_v140.txt silo
../bin/icgen_serial test_RSH/pf_rsh_v140_BC.txt silo
../bin/icgen_serial test_RSH/pf_rsh_chem_v140_BC.txt silo
#
../bin/pion_serial IC_RSH2D_n128_v140_nochem.silo   5 1 redirect=../results/msg_RSH2D_n128_v140_nochem &
../bin/pion_serial IC_RSH2D_n128_v140_nochemBC.silo 5 1 redirect=../results/msg_RSH2D_n128_v140_nochemBC &
../bin/pion_serial IC_RSH2D_n128_v140_chem.silo   5 1 redirect=../results/msg_RSH2D_n128_v140_chem &
../bin/pion_serial IC_RSH2D_n128_v140_chemBC.silo 5 1 redirect=../results/msg_RSH2D_n128_v140_chemBC &
exit

../bin/icgen_serial test_RSH/pf_rsh_v100.txt silo
../bin/icgen_serial test_RSH/pf_rsh_v120.txt silo
../bin/icgen_serial test_RSH/pf_rsh_v130.txt silo
../bin/icgen_serial test_RSH/pf_rsh_v140.txt silo
../bin/icgen_serial test_RSH/pf_rsh_v150.txt silo
../bin/icgen_serial test_RSH/pf_rsh_v170.txt silo
#
../bin/pion_serial IC_RSH2D_n128_v100_nochem.silo 5 1 redirect=../results/msg_RSH2D_n128_v100 &
../bin/pion_serial IC_RSH2D_n128_v120_nochem.silo 5 1 redirect=../results/msg_RSH2D_n128_v120 &
wait
../bin/pion_serial IC_RSH2D_n128_v130_nochem.silo 5 1 redirect=../results/msg_RSH2D_n128_v130 &
../bin/pion_serial IC_RSH2D_n128_v140_nochem.silo 5 1 redirect=../results/msg_RSH2D_n128_v140 &
wait
../bin/pion_serial IC_RSH2D_n128_v150_nochem.silo 5 1 redirect=../results/msg_RSH2D_n128_v150 &
../bin/pion_serial IC_RSH2D_n128_v170_nochem.silo 5 1 redirect=../results/msg_RSH2D_n128_v170 &
wait
#
exit

####################################
## Testing MHD code for errors!!! ##
## 2009-12                        ##
####################################
../bin/pion_serial_vanAlbadaMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp0 redirect=../results/FL_LimVA_oa2av10_Lin_ \
cfl=0.4 artvisc=0.10 solver=1 op_criterion=0 opfreq=0 &
../bin/pion_serial_vanAlbadaMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp1 redirect=../results/FL_LimVA_oa2av01_Lin_ \
cfl=0.4 artvisc=0.01 solver=1 op_criterion=0 opfreq=0 &
../bin/pion_serial_vanAlbadaMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp1 redirect=../results/FL_LimVA_oa2av00_Lin_ \
cfl=0.4 artvisc=0.00 solver=1 op_criterion=0 opfreq=0 &
wait
../bin/pion_serial_vanAlbadaMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp0 redirect=../results/FL_LimVA_oa2av10_Roe_ \
cfl=0.4 artvisc=0.10 solver=4 op_criterion=0 opfreq=0 &
../bin/pion_serial_vanAlbadaMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp1 redirect=../results/FL_LimVA_oa2av01_Roe_ \
cfl=0.4 artvisc=0.01 solver=4 op_criterion=0 opfreq=0 &
../bin/pion_serial_vanAlbadaMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp1 redirect=../results/FL_LimVA_oa2av00_Roe_ \
cfl=0.4 artvisc=0.00 solver=4 op_criterion=0 opfreq=0 &
wait
#
../bin/pion_serial_minmodMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp0 redirect=../results/FL_LimMM_oa2av10_Lin_ \
cfl=0.4 artvisc=0.10 solver=1 op_criterion=0 opfreq=0 &
../bin/pion_serial_minmodMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp1 redirect=../results/FL_LimMM_oa2av01_Lin_ \
cfl=0.4 artvisc=0.01 solver=1 op_criterion=0 opfreq=0 &
../bin/pion_serial_minmodMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp1 redirect=../results/FL_LimMM_oa2av00_Lin_ \
cfl=0.4 artvisc=0.00 solver=1 op_criterion=0 opfreq=0 &
wait
../bin/pion_serial_minmodMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp0 redirect=../results/FL_LimMM_oa2av10_Roe_ \
cfl=0.4 artvisc=0.10 solver=4 op_criterion=0 opfreq=0 &
../bin/pion_serial_minmodMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp1 redirect=../results/FL_LimMM_oa2av01_Roe_ \
cfl=0.4 artvisc=0.01 solver=4 op_criterion=0 opfreq=0 &
../bin/pion_serial_minmodMHD IC_FieldLoop200.silo 5 1 outfile=../results/FLtemp1 redirect=../results/FL_LimMM_oa2av00_Roe_ \
cfl=0.4 artvisc=0.00 solver=4 op_criterion=0 opfreq=0 &
wait

exit

make -f Makefile.serial.code; make -f Makefile.serial.icgenerator
#
../bin/pion_serial IC_FieldLoop200.silo 5 1 op_criterion=0 opfreq=200 cfl=0.4 ooa=1 artvisc=0   solver=1 outfile=mhdtest/FLfav_oa1av0
../bin/pion_serial IC_FieldLoop200.silo 5 1 op_criterion=0 opfreq=200 cfl=0.4 ooa=1 artvisc=0.1 solver=1 outfile=mhdtest/FLfav_oa1av1
../bin/pion_serial IC_FieldLoop200.silo 5 1 op_criterion=0 opfreq=200 cfl=0.4 ooa=2 artvisc=0   solver=1 outfile=mhdtest/FLfav_oa2av0
../bin/pion_serial IC_FieldLoop200.silo 5 1 op_criterion=0 opfreq=200 cfl=0.4 ooa=2 artvisc=0.1 solver=1 outfile=mhdtest/FLfav_oa2av1
#
../bin/pion_serial IC_FieldLoop200.silo 5 1 op_criterion=0 opfreq=200 cfl=0.4 ooa=1 artvisc=0   solver=4 outfile=mhdtest/FLfav_oa1av0roe
../bin/pion_serial IC_FieldLoop200.silo 5 1 op_criterion=0 opfreq=200 cfl=0.4 ooa=1 artvisc=0.1 solver=4 outfile=mhdtest/FLfav_oa1av1roe
../bin/pion_serial IC_FieldLoop200.silo 5 1 op_criterion=0 opfreq=200 cfl=0.4 ooa=2 artvisc=0   solver=4 outfile=mhdtest/FLfav_oa2av0roe
../bin/pion_serial IC_FieldLoop200.silo 5 1 op_criterion=0 opfreq=200 cfl=0.4 ooa=2 artvisc=0.1 solver=4 outfile=mhdtest/FLfav_oa2av1roe
#
../bin/pion_serial IC_OrszagTang_n256_b3.33m1.0.silo 5 1 opfreq=250 cfl=0.4 solver=1 ooa=1 artvisc=0   outfile=mhdtest/OTVfav_oa1av0
../bin/pion_serial IC_OrszagTang_n256_b3.33m1.0.silo 5 1 opfreq=250 cfl=0.4 solver=1 ooa=1 artvisc=0.1 outfile=mhdtest/OTVfav_oa1av1
../bin/pion_serial IC_OrszagTang_n256_b3.33m1.0.silo 5 1 opfreq=250 cfl=0.4 solver=1 ooa=2 artvisc=0   outfile=mhdtest/OTVfav_oa2av0
../bin/pion_serial IC_OrszagTang_n256_b3.33m1.0.silo 5 1 opfreq=250 cfl=0.4 solver=1 ooa=2 artvisc=0.1 outfile=mhdtest/OTVfav_oa2av1
#
exit


####################################
../bin/icgen_serial ../test_problems/test_ShockTubes/pf_st_falle09_hires.txt silo
../bin/pion_serial IC_Falle10kFS.silo  5 1 outfile=/mnt/local/jm/temp_sims/code_test_dir/testRoe_Falle10kFS  opfreq=1000 redirect=/mnt/local/jm/temp_sims/code_test_dir/msg_Falle10kFS_ artvisc=0.1 ooa=2 solver=4
../bin/pion_serial IC_Falle10kFS.silo  5 1 outfile=/mnt/local/jm/temp_sims/code_test_dir/testHyb_Falle10kFS  opfreq=1000 redirect=/mnt/local/jm/temp_sims/code_test_dir/msg_Falle10kFS_ artvisc=0.1 ooa=2
exit

#########################################
# 2009-11-26 testing rocket effect
#########################################
../bin/icgen_serial rocket_effect/pf_clump_d0100r05.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d0250r05.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d0500r05.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d1000r05.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d2000r05.txt silo
#
../bin/pion_serial IC_clump_dens0100_rad05.silo 5 1 outfile=../results/rocket/clump_dens0100_rad05_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0250_rad05.silo 5 1 outfile=../results/rocket/clump_dens0250_rad05_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0500_rad05.silo 5 1 outfile=../results/rocket/clump_dens0500_rad05_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens1000_rad05.silo 5 1 outfile=../results/rocket/clump_dens1000_rad05_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens2000_rad05.silo 5 1 outfile=../results/rocket/clump_dens2000_rad05_c11 solver=3 cooling=11 finishtime=9.48e12
#
../bin/pion_serial IC_clump_dens0100_rad05.silo 5 1 outfile=../results/rocket/clump_dens0100_rad05_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0250_rad05.silo 5 1 outfile=../results/rocket/clump_dens0250_rad05_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0500_rad05.silo 5 1 outfile=../results/rocket/clump_dens0500_rad05_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens1000_rad05.silo 5 1 outfile=../results/rocket/clump_dens1000_rad05_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens2000_rad05.silo 5 1 outfile=../results/rocket/clump_dens2000_rad05_c15 solver=3 cooling=15 finishtime=9.48e12
#
####
#
../bin/icgen_serial rocket_effect/pf_clump_d0025r10.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d0063r10.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d0125r10.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d0250r10.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d0500r10.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d1000r10.txt silo
../bin/icgen_serial rocket_effect/pf_clump_d2000r10.txt silo
#
../bin/pion_serial IC_clump_dens0025_rad10.silo 5 1 outfile=../results/rocket/clump_dens0025_rad10_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0063_rad10.silo 5 1 outfile=../results/rocket/clump_dens0063_rad10_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0125_rad10.silo 5 1 outfile=../results/rocket/clump_dens0125_rad10_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0250_rad10.silo 5 1 outfile=../results/rocket/clump_dens0250_rad10_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0500_rad10.silo 5 1 outfile=../results/rocket/clump_dens0500_rad10_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens1000_rad10.silo 5 1 outfile=../results/rocket/clump_dens1000_rad10_c11 solver=3 cooling=11 finishtime=9.48e12
../bin/pion_serial IC_clump_dens2000_rad10.silo 5 1 outfile=../results/rocket/clump_dens2000_rad10_c11 solver=3 cooling=11 finishtime=9.48e12
#
../bin/pion_serial IC_clump_dens0025_rad10.silo 5 1 outfile=../results/rocket/clump_dens0025_rad10_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0063_rad10.silo 5 1 outfile=../results/rocket/clump_dens0063_rad10_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0125_rad10.silo 5 1 outfile=../results/rocket/clump_dens0125_rad10_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0250_rad10.silo 5 1 outfile=../results/rocket/clump_dens0250_rad10_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens0500_rad10.silo 5 1 outfile=../results/rocket/clump_dens0500_rad10_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens1000_rad10.silo 5 1 outfile=../results/rocket/clump_dens1000_rad10_c15 solver=3 cooling=15 finishtime=9.48e12
../bin/pion_serial IC_clump_dens2000_rad10.silo 5 1 outfile=../results/rocket/clump_dens2000_rad10_c15 solver=3 cooling=15 finishtime=9.48e12
exit


#########################################
# 2009-11-17 testing artificial viscosity
#########################################
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0    solver=2 \
 outfile=test_artvisc/M10ang1_Exact_av00 redirect=test_artvisc/msg_M10ang1_Exact_av00_
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0.01 solver=2 \
 outfile=test_artvisc/M10ang1_Exact_av01 redirect=test_artvisc/msg_M10ang1_Exact_av01_
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0.10 solver=2 \
 outfile=test_artvisc/M10ang1_Exact_av10 redirect=test_artvisc/msg_M10ang1_Exact_av10_
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0.20 solver=2 \
 outfile=test_artvisc/M10ang1_Exact_av20 redirect=test_artvisc/msg_M10ang1_Exact_av20_
#
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0    solver=2 \
 outfile=test_artvisc/M25ang1_Exact_av00 redirect=test_artvisc/msg_M25ang1_Exact_av00_
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0.01 solver=2  \
 outfile=test_artvisc/M25ang1_Exact_av01 redirect=test_artvisc/msg_M25ang1_Exact_av01_
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0.10 solver=2  \
 outfile=test_artvisc/M25ang1_Exact_av10 redirect=test_artvisc/msg_M25ang1_Exact_av10_
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0.20 solver=2  \
 outfile=test_artvisc/M25ang1_Exact_av20 redirect=test_artvisc/msg_M25ang1_Exact_av20_
exit
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0    solver=4 \
 outfile=test_artvisc/M10ang1_Roe_av00 redirect=test_artvisc/msg_M10ang1_Roe_av00_
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0.01 solver=4 \
 outfile=test_artvisc/M10ang1_Roe_av01 redirect=test_artvisc/msg_M10ang1_Roe_av01_
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0.10 solver=4 \
 outfile=test_artvisc/M10ang1_Roe_av10 redirect=test_artvisc/msg_M10ang1_Roe_av10_
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0.20 solver=4 \
 outfile=test_artvisc/M10ang1_Roe_av20 redirect=test_artvisc/msg_M10ang1_Roe_av20_
#
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0    solver=4 \
 outfile=test_artvisc/M25ang1_Roe_av00 redirect=test_artvisc/msg_M25ang1_Roe_av00_
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0.01 solver=4  \
 outfile=test_artvisc/M25ang1_Roe_av01 redirect=test_artvisc/msg_M25ang1_Roe_av01_
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0.10 solver=4  \
 outfile=test_artvisc/M25ang1_Roe_av10 redirect=test_artvisc/msg_M25ang1_Roe_av10_
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0.20 solver=4  \
 outfile=test_artvisc/M25ang1_Roe_av20 redirect=test_artvisc/msg_M25ang1_Roe_av20_
exit
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0 \
 outfile=test_artvisc/M10ang1_av00 redirect=test_artvisc/msg_M10ang1_av00_
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0.01 \
 outfile=test_artvisc/M10ang1_av01 redirect=test_artvisc/msg_M10ang1_av01_
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0.10 \
 outfile=test_artvisc/M10ang1_av10 redirect=test_artvisc/msg_M10ang1_av10_
../bin/pion_serial test_artvisc/IC_M10shock_angle01deg.silo 5 1 artvisc=0.20 \
 outfile=test_artvisc/M10ang1_av20 redirect=test_artvisc/msg_M10ang1_av20_
#
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0    \
 outfile=test_artvisc/M25ang1_av00 redirect=test_artvisc/msg_M25ang1_av00_
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0.01  \
 outfile=test_artvisc/M25ang1_av01 redirect=test_artvisc/msg_M25ang1_av01_
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0.10  \
 outfile=test_artvisc/M25ang1_av10 redirect=test_artvisc/msg_M25ang1_av10_
../bin/pion_serial test_artvisc/IC_M25shock_angle01deg.silo 5 1 artvisc=0.20  \
 outfile=test_artvisc/M25ang1_av20 redirect=test_artvisc/msg_M25ang1_av20_
exit

../bin/pion_serial test_artvisc/IC_M25shock_ang01deg.silo 5 1 artvisc=0    finishtime=3.16e10 opfreq_time=100.0 op_criterion=1 \
 outfile=test_artvisc/M25ang1_av00 redirect=test_artvisc/msg_M25ang1_av00_
../bin/pion_serial test_artvisc/IC_M25shock_ang01deg.silo 5 1 artvisc=0.01 finishtime=3.16e10 opfreq_time=100.0 op_criterion=1 \
 outfile=test_artvisc/M25ang1_av01 redirect=test_artvisc/msg_M25ang1_av01_
../bin/pion_serial test_artvisc/IC_M25shock_ang01deg.silo 5 1 artvisc=0.10 finishtime=3.16e10 opfreq_time=100.0 op_criterion=1 \
 outfile=test_artvisc/M25ang1_av10 redirect=test_artvisc/msg_M25ang1_av10_
#
../bin/pion_serial test_artvisc/IC_M25shock_ang01deg.silo 5 1 artvisc=0    outfile=test_artvisc/M25ang1_Roe_av00 \
finishtime=3.16e10 opfreq_time=100.0 op_criterion=1 solver=4 redirect=test_artvisc/msg_M25ang1_Roe_av00_
../bin/pion_serial test_artvisc/IC_M25shock_ang01deg.silo 5 1 artvisc=0.01 outfile=test_artvisc/M25ang1_Roe_av01 \
finishtime=3.16e10 opfreq_time=100.0 op_criterion=1 solver=4 redirect=test_artvisc/msg_M25ang1_Roe_av01_
../bin/pion_serial test_artvisc/IC_M25shock_ang01deg.silo 5 1 artvisc=0.10 outfile=test_artvisc/M25ang1_Roe_av10 \
finishtime=3.16e10 opfreq_time=100.0 op_criterion=1 solver=4 redirect=test_artvisc/msg_M25ang1_Roe_av10_
exit




#################
# 2009-10-21/22 testing Roe MHD flux solver:
#################

../bin/pion_serial pf_st_falle07.txt 1 1 solver=1 outfile=../results/FALLEav07     finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2 >tmp.txt
../bin/pion_serial pf_st_falle07.txt 1 1 solver=4 outfile=../results/FALLEav07_Roe finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2 >>tmp.txt
../bin/pion_serial pf_st_falle09.txt 1 1 solver=1 outfile=../results/FALLEav09     finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2 >>tmp.txt
../bin/pion_serial pf_st_falle09.txt 1 1 solver=4 outfile=../results/FALLEav09_Roe finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2  >>tmp.txt
../bin/pion_serial pf_st_falle10.txt 1 1 solver=1 outfile=../results/FALLEav10     finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2  >>tmp.txt
../bin/pion_serial pf_st_falle10.txt 1 1 solver=4 outfile=../results/FALLEav10_Roe finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2  >>tmp.txt
../bin/pion_serial pf_st_falle11.txt 1 1 solver=1 outfile=../results/FALLEav11     finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2  >>tmp.txt
../bin/pion_serial pf_st_falle11.txt 1 1 solver=4 outfile=../results/FALLEav11_Roe finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2  >>tmp.txt
../bin/pion_serial pf_st_falle12.txt 1 1 solver=1 outfile=../results/FALLEav12     finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2  >>tmp.txt
../bin/pion_serial pf_st_falle12.txt 1 1 solver=4 outfile=../results/FALLEav12_Roe finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2  >>tmp.txt
../bin/pion_serial pf_st_falle13.txt 1 1 solver=1 outfile=../results/FALLEav13     finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2  >>tmp.txt
../bin/pion_serial pf_st_falle13.txt 1 1 solver=4 outfile=../results/FALLEav13_Roe finishtime=0.2 op_criterion=1 opfreq_time=1.582278481e-9 ooa=2 artvisc=0.15 eqntype=2  >>tmp.txt
exit



#################
# 2009-10-21/22 testing new solver structure:
#################
valgrind  ../bin/pion_serial pf_st_toro1.txt 1 1 solver=3 ooa=2 artvisc=0.1 op_criterion=0 outfile=../results/TT1_RS3_OA2_av1
valgrind  ../bin/pion_serial pf_st_toro2.txt 1 1 solver=3 ooa=2 artvisc=0.1 op_criterion=0 outfile=../results/TT2_RS3_OA2_av1
valgrind  ../bin/pion_serial pf_st_toro3.txt 1 1 solver=3 ooa=2 artvisc=0.1 op_criterion=0 outfile=../results/TT3_RS3_OA2_av1
valgrind  ../bin/pion_serial pf_st_toro4.txt 1 1 solver=3 ooa=2 artvisc=0.1 op_criterion=0 outfile=../results/TT4_RS3_OA2_av1
valgrind  ../bin/pion_serial pf_st_toro5.txt 1 1 solver=3 ooa=2 artvisc=0.1 op_criterion=0 outfile=../results/TT5_RS3_OA2_av1
exit


make -f Makefile.serial.icgenerator; make -f Makefile.serial.code
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RSlinear_OA1_av1 op_criterion=0 opfreq=100 solver=1 ooa=1 artvisc=0.1 eqntype=1 redirect=msg_dmr_RSlinear_OA1_av1_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RSlinear_OA2_av1 op_criterion=0 opfreq=100 solver=1 ooa=2 artvisc=0.1 eqntype=1 redirect=msg_dmr_RSlinear_OA2_av1_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RSlinear_OA1_av0 op_criterion=0 opfreq=100 solver=1 ooa=1 artvisc=0.0 eqntype=1 redirect=msg_dmr_RSlinear_OA1_av0_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RSlinear_OA2_av0 op_criterion=0 opfreq=100 solver=1 ooa=2 artvisc=0.0 eqntype=1 redirect=msg_dmr_RSlinear_OA2_av0_
#
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RSexact_OA1_av1 op_criterion=0 opfreq=100 solver=2 ooa=1 artvisc=0.1 eqntype=1 redirect=msg_dmr_RSexact_OA1_av1_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RSexact_OA2_av1 op_criterion=0 opfreq=100 solver=2 ooa=2 artvisc=0.1 eqntype=1 redirect=msg_dmr_RSexact_OA2_av1_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RSexact_OA1_av0 op_criterion=0 opfreq=100 solver=2 ooa=1 artvisc=0.0 eqntype=1 redirect=msg_dmr_RSexact_OA1_av0_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RSexact_OA2_av0 op_criterion=0 opfreq=100 solver=2 ooa=2 artvisc=0.0 eqntype=1 redirect=msg_dmr_RSexact_OA2_av0_
#
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RShybrid_OA1_av1 op_criterion=0 opfreq=100 solver=3 ooa=1 artvisc=0.1 eqntype=1 redirect=msg_dmr_RShybrid_OA1_av1_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RShybrid_OA2_av1 op_criterion=0 opfreq=100 solver=3 ooa=2 artvisc=0.1 eqntype=1 redirect=msg_dmr_RShybrid_OA2_av1_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RShybrid_OA1_av0 op_criterion=0 opfreq=100 solver=3 ooa=1 artvisc=0.0 eqntype=1 redirect=msg_dmr_RShybrid_OA1_av0_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_RShybrid_OA2_av0 op_criterion=0 opfreq=100 solver=3 ooa=2 artvisc=0.0 eqntype=1 redirect=msg_dmr_RShybrid_OA2_av0_
#
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_LF_OA1_av0 op_criterion=0 opfreq=100 solver=0 ooa=1 artvisc=0.0 eqntype=1 redirect=msg_dmr_LF_OA1_av0_
../bin/pion_serial IC_DMRm10t60_n512.silo 5 1 outfile=../results/dmr_new_LF_OA2_av0 op_criterion=0 opfreq=100 solver=0 ooa=2 artvisc=0.0 eqntype=1 redirect=msg_dmr_LF_OA2_av0_
exit


###############################################
# 31/8/09 Testing the Even Better Henney Cooling
# Model, to see if I can avoid crashes...
###############################################
make -f Makefile.serial.code

../bin/pion_serial IC_S80S2D_n256_Euler.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_Euler_COOLHEATav2_  \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_Euler_COOLHEATav2  artvisc=0.2 &
../bin/pion_serial IC_S80S2D_n256_MHDbeta1p5.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_MHDbeta1p5_COOLHEATav2_ \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_MHDbeta1p5_COOLHEATav2  artvisc=0.2 &
../bin/pion_serial IC_S80S2D_n256_MHDbeta0p15.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_MHDbeta0p15_COOLHEATav2_ \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_MHDbeta0p15_COOLHEATav2  artvisc=0.2 &
../bin/pion_serial IC_S80S2D_n256_MHDbeta0p015.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_MHDbeta0p015_COOLHEATav2_ \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_MHDbeta0p015_COOLHEATav2  artvisc=0.2 &
wait

../bin/pion_serial IC_S80S2D_n256_Euler.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_Euler_COOLHEATav1_ \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_Euler_COOLHEATav1  artvisc=0.1 &
../bin/pion_serial IC_S80S2D_n256_MHDbeta1p5.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_MHDbeta1p5_COOLHEATav1_ \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_MHDbeta1p5_COOLHEATav1  artvisc=0.1 &
../bin/pion_serial IC_S80S2D_n256_MHDbeta0p15.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_MHDbeta0p15_COOLHEATav1_ \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_MHDbeta0p15_COOLHEATav1  artvisc=0.1 &
../bin/pion_serial IC_S80S2D_n256_MHDbeta0p015.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_MHDbeta0p015_COOLHEATav1_ \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_MHDbeta0p015_COOLHEATav1  artvisc=0.1 &
wait

../bin/pion_serial IC_S80S2D_n512_Euler.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n512_Euler_COOLHEATav2_  \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n512_Euler_COOLHEATav2  artvisc=0.2 &
../bin/pion_serial IC_S80S2D_n1024_Euler.silo 5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n1024_Euler_COOLHEATav2_ \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n1024_Euler_COOLHEATav2  artvisc=0.2 &
wait

../bin/pion_serial IC_S80S2D_n512_Euler.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n512_Euler_COOLHEATav1_  \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n512_Euler_COOLHEATav1   artvisc=0.1 &
../bin/pion_serial IC_S80S2D_n1024_Euler.silo 5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n1024_Euler_COOLHEATav1_ \
outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n1024_Euler_COOLHEATav1  artvisc=0.1 &
wait
exit

###############################################
# 29/8/09 Test run for modified code structure.
###############################################
#../bin/pion_serial_orig IC_S80S2D_n256_Euler.silo  5 1 op_criterion=0 opfreq=50 finishtime=8.3e10 outfile=CODETEST_EULER_OLD &
../bin/pion_serial      IC_S80S2D_n256_Euler.silo  5 1 op_criterion=0 opfreq=50 finishtime=8.3e10 outfile=CODETEST_EULER_NEW &
#../bin/pion_serial_orig IC_S80S2D_n256_MHDbeta1p5.silo  5 1 op_criterion=0 opfreq=50 finishtime=8.3e10 outfile=CODETEST_MHDbeta1p5_OLD &
../bin/pion_serial      IC_S80S2D_n256_MHDbeta1p5.silo  5 1 op_criterion=0 opfreq=50 finishtime=8.3e10 outfile=CODETEST_MHDbeta1p5_NEW &
wait
cd /mnt/local/jm/local_libs/silo_gcc_hdf5/bin/
./silodiff ~/active/projects/uniform_grid_code/trunk/bin_serial/CODETEST_EULER_NEW.000050.silo \
           ~/active/projects/uniform_grid_code/trunk/bin_serial/CODETEST_EULER_OLD.000050.silo 
./silodiff ~/active/projects/uniform_grid_code/trunk/bin_serial/CODETEST_MHDbeta1p5_NEW.000043.silo \
           ~/active/projects/uniform_grid_code/trunk/bin_serial/CODETEST_MHDbeta1p5_OLD.000043.silo 
cd /home/jmackey/active/projects/uniform_grid_code/trunk/bin_serial/
exit

###############################################
# 26/8/09 Testing the Better Henney Cooling
# Model, to see if I can avoid crashes...
###############################################
../bin/pion_serial IC_S80S2D_n256_Euler.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_Euler_  outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_Euler &
../bin/pion_serial IC_S80S2D_n512_Euler.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n512_Euler_  outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n512_Euler &
../bin/pion_serial IC_S80S2D_n1024_Euler.silo 5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n1024_Euler_ outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n1024_Euler &
wait
../bin/pion_serial_heat26 IC_S80S2D_n256_Euler.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n256_Euler_heat26_  outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n256_Euler_heat26 &
../bin/pion_serial_heat26 IC_S80S2D_n512_Euler.silo  5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n512_Euler_heat26_  outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n512_Euler_heat26 &
../bin/pion_serial_heat26 IC_S80S2D_n1024_Euler.silo 5 1 redirect=HenneyCoolingTesting/msg/msg_S80S2D_n1024_Euler_heat26_ outfile=/mnt/local/jm/mysims/HenneyCoolingTests/S80S2D_n1024_Euler_heat26 &
wait
exit


###############################################
# 30/7/09 More testing 2D MHD with Henney's cooling.
###############################################
../bin/pion_serial IC_TEST3_n128_Henneyetal09.silo 5 1    opfreq=050 redirect=../results/test_cool_20090730/msg_TEST3_n128eul_c16 outfile=../results/test_cool_20090730/TEST3_n128eul_c16
../bin/pion_serial IC_TEST3_n128mhd_Henneyetal09.silo 5 1 opfreq=050 redirect=../results/test_cool_20090730/msg_TEST3_n128mhd_c16 outfile=../results/test_cool_20090730/TEST3_n128mhd_c16
../bin/pion_serial IC_TEST3_n256_Henneyetal09.silo 5 1    opfreq=100 redirect=../results/test_cool_20090730/msg_TEST3_n256eul_c16 outfile=../results/test_cool_20090730/TEST3_n256eul_c16
../bin/pion_serial IC_TEST3_n256mhd_Henneyetal09.silo 5 1 opfreq=100 redirect=../results/test_cool_20090730/msg_TEST3_n256mhd_c16 outfile=../results/test_cool_20090730/TEST3_n256mhd_c16
../bin/pion_serial IC_TEST3_n512_Henneyetal09.silo 5 1    opfreq=200 redirect=../results/test_cool_20090730/msg_TEST3_n512eul_c16 outfile=../results/test_cool_20090730/TEST3_n512eul_c16
../bin/pion_serial IC_TEST3_n512mhd_Henneyetal09.silo 5 1 opfreq=200 redirect=../results/test_cool_20090730/msg_TEST3_n512mhd_c16 outfile=../results/test_cool_20090730/TEST3_n512mhd_c16
../bin/pion_serial IC_TEST3_n1024_Henneyetal09.silo 5 1    opfreq=400 redirect=../results/test_cool_20090730/msg_TEST3_n1024eul_c16 outfile=../results/test_cool_20090730/TEST3_n1024eul_c16
../bin/pion_serial IC_TEST3_n1024mhd_Henneyetal09.silo 5 1 opfreq=400 redirect=../results/test_cool_20090730/msg_TEST3_n1024mhd_c16 outfile=../results/test_cool_20090730/TEST3_n1024mhd_c16
../bin/pion_serial IC_TEST3_n2048_Henneyetal09.silo 5 1    opfreq=800 redirect=../results/test_cool_20090730/msg_TEST3_n2048eul_c16 outfile=../results/test_cool_20090730/TEST3_n2048eul_c16
../bin/pion_serial IC_TEST3_n2048mhd_Henneyetal09.silo 5 1 opfreq=800 redirect=../results/test_cool_20090730/msg_TEST3_n2048mhd_c16 outfile=../results/test_cool_20090730/TEST3_n2048mhd_c16


###############################################
# 27/7/09 Testing 2D MHD with Henney cooling. #
###############################################
make -f Makefile.serial.code
#../bin/pion_serial IC_R2d_n256_Henneyetal09.silo 5 1 outfile=../results/T2D256_cfl10_Henneyetal09_C15 cooling=15 cfl=0.10 redirect=msg_Henney_n256_cfl10_C15_ opfreq_time=5000.0 &
#../bin/pion_serial IC_R2d_n512_Henneyetal09.silo 5 1 outfile=../results/T2D512_cfl30_Henneyetal09_C15 cooling=15 cfl=0.30 redirect=msg_Henney_n512_cfl30_C15_ opfreq_time=5000.0 &
wait
../bin/pion_serial IC_R2d_n256_Henneyetal09.silo 5 1 outfile=../results/T2Dv2_n256_cfl05_Henneyetal09 cfl=0.05 redirect=msg_HenneyV2_n256_cfl05_ opfreq_time=5000.0 &
../bin/pion_serial IC_R2d_n512_Henneyetal09.silo 5 1 outfile=../results/T2Dv2_n512_cfl05_Henneyetal09 cfl=0.05 redirect=msg_HenneyV2_n512_cfl05_ opfreq_time=5000.0 &
wait

../bin/pion_serial IC_R2d_n1024_Henneyetal09.silo 5 1 outfile=../results/T2Dv2_n1024_cfl10_Henneyetal09 cfl=0.10 redirect=msg_HenneyV2_n1024_cfl10_ opfreq_time=5000.0 &
../bin/pion_serial IC_R2d_n1024_Henneyetal09.silo 5 1 outfile=../results/T2Dv2_n1024_cfl05_Henneyetal09 cfl=0.05 redirect=msg_HenneyV2_n1024_cfl05_ opfreq_time=5000.0 &

wait
#
#../bin/pion_serial IC_R2d_n1024_Henneyetal09.silo 5 1 outfile=../results/T2D1024_cfl05_Henneyetal09 cfl=0.05 redirect=msg_Henney_n1024_cfl05_ opfreq_time=5000.0

exit





#################################
# 16/5/09 2D SS instabilities   #
#################################
../bin/pion_serial IC_ss2d_ref_n512_noisy.silo 5 1 outfile=../results/r2dtesting/ss2d_ref_n512_noisy_c11 \
cooling=11 finishtime=9.48e12 redirect=../results/r2dtesting/msg_ss2d_ref_n512_noisy_c11 \
op_criterion=1 opfreq_time=10000.0 &
../bin/pion_serial IC_ss2d_ref_n512_noisy.silo 5 1 outfile=../results/r2dtesting/ss2d_ref_n512_noisy_c12 \
cooling=12 finishtime=9.48e12 redirect=../results/r2dtesting/msg_ss2d_ref_n512_noisy_c12  \
op_criterion=1 opfreq_time=10000.0 &
wait

#../bin/pion_serial IC_ss2d_out_n512_noisy.silo 5 1 outfile=../results/r2dtesting/ss2d_out_n512_noisy_c11 \
#cooling=11 finishtime=9.48e12 redirect=../results/r2dtesting/msg_ss2d_out_n512_noisy_c11 &
#../bin/pion_serial IC_ss2d_out_n512_noisy.silo 5 1 outfile=../results/r2dtesting/ss2d_out_n512_noisy_c12 \
#cooling=12 finishtime=9.48e12 redirect=../results/r2dtesting/msg_ss2d_out_n512_noisy_c12 &
#wait

exit

#################################
# 6/5/09 testing clump positions#
#################################
make -f Makefile.serial.code
#../bin/pion_serial IC_R2dLine0.silo 5 1 outfile=Line_T0
#../bin/pion_serial IC_R2dLine10.silo 5 1 outfile=Line_T10
#../bin/pion_serial IC_R2dLine20.silo 5 1 outfile=Line_T20
#../bin/pion_serial IC_R2dLine30.silo 5 1 outfile=Line_T30
../bin/pion_serial IC_R2dLine_theta0.silo  5 1 outfile=Line_clumps2d/Line_T00_cool cooling=12 redirect=Line_clumps2d/msgT00_cool_ &
../bin/pion_serial IC_R2dLine_theta10.silo 5 1 outfile=Line_clumps2d/Line_T10_cool cooling=12 redirect=Line_clumps2d/msgT10_cool_ &
wait
../bin/pion_serial IC_R2dLine_theta20.silo 5 1 outfile=Line_clumps2d/Line_T20_cool cooling=12 redirect=Line_clumps2d/msgT20_cool_ &
../bin/pion_serial IC_R2dLine_theta30.silo 5 1 outfile=Line_clumps2d/Line_T30_cool cooling=12 redirect=Line_clumps2d/msgT30_cool_ &
wait
exit


###########################################
# testing code. 13/2/09
###########################################
make -f Makefile.serial.code
#
# short sims -- just a few steps
#
../bin/pion_serial ../results/diff_tests_200902/IC_random3d_n88_serial.silo 5 1 \
outfile=../results/diff_tests_200902/random3d_epona_serial \
redirect=../results/diff_tests_200902/msg_random3d_epona_serial_ \
op_criterion=0 opfreq=1 finishtime=1.30e10
#
../bin/pion_serial ../results/diff_tests_200902/IC_random2d_n240_serial.silo 5 1 \
outfile=../results/diff_tests_200902/random2d_epona_serial \
redirect=../results/diff_tests_200902/msg_random2d_epona_serial_ \
op_criterion=0 opfreq=1 finishtime=1.30e10
#
## longer sims
#
../bin/pion_serial ../results/diff_tests_200902/IC_random3d_n88_serial.silo 5 1 \
outfile=../results/diff_tests_200902/random3d_epona_serial \
redirect=../results/diff_tests_200902/msg_random3d_epona_serial_ opfreq_time=5000.0
#
../bin/pion_serial ../results/diff_tests_200902/IC_random2d_n240_serial.silo 5 1 \
outfile=../results/diff_tests_200902/random2d_epona_serial \
redirect=../results/diff_tests_200902/msg_random2d_epona_serial_ opfreq_time=5000.0
#
exit
###########################################


