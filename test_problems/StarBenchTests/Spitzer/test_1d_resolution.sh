#!/bin/bash



../../../icgen_serial params_Spitzer1D_n0080.txt silo
../../../icgen_serial params_Spitzer1D_n0160.txt silo
../../../icgen_serial params_Spitzer1D_n0320.txt silo
../../../icgen_serial params_Spitzer1D_n0640.txt silo
../../../icgen_serial params_Spitzer1D_n1280.txt silo

mkdir data
mkdir data_ascii

../../../pion_serial IC_Spitzer1D_n0080.silo 5 1 optype=silo outfile=data/S1D_n0080 redirect=log_S1Dn0080 # &
#../../../pion_serial IC_Spitzer1D_n0080.silo 5 1 optype=text outfile=data_ascii/S1D_n0080 redirect=log_junk &
#wait
../../../pion_serial IC_Spitzer1D_n0160.silo 5 1 optype=silo outfile=data/S1D_n0160 redirect=log_S1Dn0160 #&
#../../../pion_serial IC_Spitzer1D_n0160.silo 5 1 optype=text outfile=data_ascii/S1D_n0160 redirect=log_junk &
#wait
../../../pion_serial IC_Spitzer1D_n0320.silo 5 1 optype=silo outfile=data/S1D_n0320 redirect=log_S1Dn0320 #&
#../../../pion_serial IC_Spitzer1D_n0320.silo 5 1 optype=text outfile=data_ascii/S1D_n0320 redirect=log_junk &
#wait
../../../pion_serial IC_Spitzer1D_n0640.silo 5 1 optype=silo outfile=data/S1D_n0640 redirect=log_S1Dn0640 #&
#../../../pion_serial IC_Spitzer1D_n0640.silo 5 1 optype=text outfile=data_ascii/S1D_n0640 redirect=log_junk &
#wait
../../../pion_serial IC_Spitzer1D_n1280.silo 5 1 optype=silo outfile=data/S1D_n1280 redirect=log_S1Dn1280 #&
#../../../pion_serial IC_Spitzer1D_n1280.silo 5 1 optype=text outfile=data_ascii/S1D_n1280 redirect=log_junk &
#wait

cd analysis
make

./Radii ../data S1D_n0080.0 . S1Dn0080
./Radii ../data S1D_n0160.0 . S1Dn0160
./Radii ../data S1D_n0320.0 . S1Dn0320
./Radii ../data S1D_n0640.0 . S1Dn0640
./Radii ../data S1D_n1280.0 . S1Dn1280


