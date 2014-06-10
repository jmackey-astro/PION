#!/bin/bash
#
# 2013.07.18 JM: Script to run 1D resolution test for HII region D-type expansion.
# 2013.07.19 JM: Modified test so that I can resolve the ionisation length.
#

../../../icgen_serial params_HII1Dv2_n00128.txt silo
../../../icgen_serial params_HII1Dv2_n00256.txt silo
../../../icgen_serial params_HII1Dv2_n00512.txt silo
../../../icgen_serial params_HII1Dv2_n01024.txt silo
../../../icgen_serial params_HII1Dv2_n02048.txt silo
../../../icgen_serial params_HII1Dv2_n04096.txt silo
../../../icgen_serial params_HII1Dv2_n08192.txt silo
../../../icgen_serial params_HII1Dv2_n16384.txt silo
../../../icgen_serial params_HII1Dv2_n32768.txt silo

data_ascii=/vol/klaipeda3/scratch/jmackey/HII1Dv2
data_silo=${data_ascii}/silo

#### TESTING ###
#### TESTING ###


mkdir ${data_ascii}
mkdir ${data_silo}
rm ${data_silo}/*.silo
rm ${data_ascii}/*.txt

../../../pion_serial IC_HII1Dv2_n00128.silo 5 1 optype=6 outfile=${data_ascii}/HII1Dv2_n00128 redirect=log00128 &
../../../pion_serial IC_HII1Dv2_n00256.silo 5 1 optype=6 outfile=${data_ascii}/HII1Dv2_n00256 redirect=log00256 &
../../../pion_serial IC_HII1Dv2_n00512.silo 5 1 optype=6 outfile=${data_ascii}/HII1Dv2_n00512 redirect=log00512 &
../../../pion_serial IC_HII1Dv2_n01024.silo 5 1 optype=6 outfile=${data_ascii}/HII1Dv2_n01024 redirect=log01024 &
../../../pion_serial IC_HII1Dv2_n02048.silo 5 1 optype=6 outfile=${data_ascii}/HII1Dv2_n02048 redirect=log02048 &
../../../pion_serial IC_HII1Dv2_n04096.silo 5 1 optype=6 outfile=${data_ascii}/HII1Dv2_n04096 redirect=log04096 &
../../../pion_serial IC_HII1Dv2_n08192.silo 5 1 optype=6 outfile=${data_ascii}/HII1Dv2_n08192 redirect=log08192 &
../../../pion_serial IC_HII1Dv2_n16384.silo 5 1 optype=6 outfile=${data_ascii}/HII1Dv2_n16384 redirect=log16384 &
../../../pion_serial IC_HII1Dv2_n32768.silo 5 1 optype=6 outfile=${data_ascii}/HII1Dv2_n32768 redirect=log32768 &
wait

mv ${data_ascii}/*.silo ${data_silo}/
cd analysis
make -j8

./Radii ${data_silo} HII1Dv2_n00128.0 . S1Dn00128
./Radii ${data_silo} HII1Dv2_n00256.0 . S1Dn00256
./Radii ${data_silo} HII1Dv2_n00512.0 . S1Dn00512
./Radii ${data_silo} HII1Dv2_n01024.0 . S1Dn01024
./Radii ${data_silo} HII1Dv2_n02048.0 . S1Dn02048
./Radii ${data_silo} HII1Dv2_n04096.0 . S1Dn04096
./Radii ${data_silo} HII1Dv2_n08192.0 . S1Dn08192
./Radii ${data_silo} HII1Dv2_n16384.0 . S1Dn16384
./Radii ${data_silo} HII1Dv2_n32768.0 . S1Dn32768


exit

