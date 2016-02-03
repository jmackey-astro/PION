#!/bin/bash
#
# 2013.07.18 JM: Script to run 1D resolution test for HII region D-type expansion.
#

../../../icgen_serial params_Spitzer1D_c2R10_n0080.txt silo
../../../icgen_serial params_Spitzer1D_c2R10_n0160.txt silo
../../../icgen_serial params_Spitzer1D_c2R10_n0320.txt silo
../../../icgen_serial params_Spitzer1D_c2R10_n0640.txt silo
../../../icgen_serial params_Spitzer1D_c2R10_n1280.txt silo
../../../icgen_serial params_Spitzer1D_c2R10_n2560.txt silo
../../../icgen_serial params_Spitzer1D_c2R10_n5120.txt silo

data_silo=/vol/klaipeda3/scratch/jmackey/Sc2R10_silo
data_ascii=/vol/klaipeda3/scratch/jmackey/Sc2R10_ascii
data_silo=./data/Sc2R10_silo
data_ascii=./data/Sc2R10_ascii

#### TESTING ###
#### TESTING ###


mkdir ${data_silo}
mkdir ${data_ascii}
rm ${data_silo}/*.silo
rm ${data_ascii}/*.txt

#../../../pion_serial IC_Spitzer1D_c2R10_n0080.silo 5 1 optype=silo outfile=${data_silo}/S1D_n0080 redirect=log_S1Dn0080  &
../../../pion_serial IC_Spitzer1D_c2R10_n0080.silo 5 1 optype=text outfile=${data_ascii}/S1D_n0080 redirect=log_junk &
#../../../pion_serial IC_Spitzer1D_c2R10_n0160.silo 5 1 optype=silo outfile=${data_silo}/S1D_n0160 redirect=log_S1Dn0160 &
../../../pion_serial IC_Spitzer1D_c2R10_n0160.silo 5 1 optype=text outfile=${data_ascii}/S1D_n0160 redirect=log_junk &
#../../../pion_serial IC_Spitzer1D_c2R10_n0320.silo 5 1 optype=silo outfile=${data_silo}/S1D_n0320 redirect=log_S1Dn0320 &
../../../pion_serial IC_Spitzer1D_c2R10_n0320.silo 5 1 optype=text outfile=${data_ascii}/S1D_n0320 redirect=log_junk &
#../../../pion_serial IC_Spitzer1D_c2R10_n0640.silo 5 1 optype=silo outfile=${data_silo}/S1D_n0640 redirect=log_S1Dn0640 &
../../../pion_serial IC_Spitzer1D_c2R10_n0640.silo 5 1 optype=text outfile=${data_ascii}/S1D_n0640 redirect=log_junk &
#../../../pion_serial IC_Spitzer1D_c2R10_n1280.silo 5 1 optype=silo outfile=${data_silo}/S1D_n1280 redirect=log_S1Dn1280 &
../../../pion_serial IC_Spitzer1D_c2R10_n1280.silo 5 1 optype=text outfile=${data_ascii}/S1D_n1280 redirect=log_junk &
#../../../pion_serial IC_Spitzer1D_c2R10_n2560.silo 5 1 optype=silo outfile=${data_silo}/S1D_n2560 redirect=log_S1Dn2560 &
../../../pion_serial IC_Spitzer1D_c2R10_n2560.silo 5 1 optype=text outfile=${data_ascii}/S1D_n2560 redirect=log_junk &
#../../../pion_serial IC_Spitzer1D_c2R10_n5120.silo 5 1 optype=silo outfile=${data_silo}/S1D_n5120 redirect=log_S1Dn5120 &
../../../pion_serial IC_Spitzer1D_c2R10_n5120.silo 5 1 optype=text outfile=${data_ascii}/S1D_n5120 redirect=log_junk &
wait

exit

cd analysis
make -j8

./Radii ${data_silo} S1D_n0080.0 . S1Dn0080
./Radii ${data_silo} S1D_n0160.0 . S1Dn0160
./Radii ${data_silo} S1D_n0320.0 . S1Dn0320
./Radii ${data_silo} S1D_n0640.0 . S1Dn0640
./Radii ${data_silo} S1D_n1280.0 . S1Dn1280
./Radii ${data_silo} S1D_n2560.0 . S1Dn2560
./Radii ${data_silo} S1D_n5120.0 . S1Dn5120


