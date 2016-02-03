#!/bin/bash
#
# 2013.06.25 JM: Script to run 1D resolution test for HII region D-type expansion.
#

data_silo=/vol/klaipeda3/scratch/jmackey/Spitzer_siloV2
data_ascii=/vol/klaipeda3/scratch/jmackey/Spitzer_asciiV2

#### TESTING ###
cd analysis
sed -i -e "s/^\/\/\#define EXCLUDE_HD_MODULE/\#define EXCLUDE_HD_MODULE/" ../../../../source/defines/functionality_flags.h
make -j8

./Radii ${data_silo} S1D_n00080.0 . S1Dn00080
./Radii ${data_silo} S1D_n00160.0 . S1Dn00160
./Radii ${data_silo} S1D_n00320.0 . S1Dn00320
./Radii ${data_silo} S1D_n00640.0 . S1Dn00640
./Radii ${data_silo} S1D_n01280.0 . S1Dn01280
./Radii ${data_silo} S1D_n02560.0 . S1Dn02560
./Radii ${data_silo} S1D_n05120.0 . S1Dn05120
./Radii ${data_silo} S1D_n10240.0 . S1Dn10240
./Radii ${data_silo} S1D_n20480.0 . S1Dn20480
./Radii ${data_silo} S1D_n40960.0 . S1Dn40960

bash plot.sh

cd -
exit
#### TESTING ###


mkdir ${data_silo}
mkdir ${data_ascii}
#rm ${data_silo}/*.silo ${data_silo}/*.txt
#rm ${data_ascii}/*.txt

cp ../../../icgen_serial .
cp ../../../pion_serial  .

./icgen_serial params_Spitzer1D_n0080.txt silo
./icgen_serial params_Spitzer1D_n0160.txt silo
./icgen_serial params_Spitzer1D_n0320.txt silo
./icgen_serial params_Spitzer1D_n0640.txt silo
./icgen_serial params_Spitzer1D_n1280.txt silo
./icgen_serial params_Spitzer1D_n2560.txt silo
./icgen_serial params_Spitzer1D_n5120.txt silo
./icgen_serial params_Spitzer1D_n10240.txt silo
./icgen_serial params_Spitzer1D_n20480.txt silo
./icgen_serial params_Spitzer1D_n40960.txt silo



./pion_serial IC_Spitzer1D_n00080.silo 5 1 optype=6 outfile=${data_silo}/S1D_n00080 redirect=log_S1Dn00080 &
./pion_serial IC_Spitzer1D_n00160.silo 5 1 optype=6 outfile=${data_silo}/S1D_n00160 redirect=log_S1Dn00160 &
./pion_serial IC_Spitzer1D_n00320.silo 5 1 optype=6 outfile=${data_silo}/S1D_n00320 redirect=log_S1Dn00320 &
./pion_serial IC_Spitzer1D_n00640.silo 5 1 optype=6 outfile=${data_silo}/S1D_n00640 redirect=log_S1Dn00640 &
./pion_serial IC_Spitzer1D_n01280.silo 5 1 optype=6 outfile=${data_silo}/S1D_n01280 redirect=log_S1Dn01280 &
./pion_serial IC_Spitzer1D_n02560.silo 5 1 optype=6 outfile=${data_silo}/S1D_n02560 redirect=log_S1Dn02560 &
./pion_serial IC_Spitzer1D_n05120.silo 5 1 optype=6 outfile=${data_silo}/S1D_n05120 redirect=log_S1Dn05120 &
./pion_serial IC_Spitzer1D_n10240.silo 5 1 optype=6 outfile=${data_silo}/S1D_n10240 redirect=log_S1Dn10240 &
./pion_serial IC_Spitzer1D_n20480.silo 5 1 optype=6 outfile=${data_silo}/S1D_n20480 redirect=log_S1Dn20480 &
./pion_serial IC_Spitzer1D_n40960.silo 5 1 optype=6 outfile=${data_silo}/S1D_n40960 redirect=log_S1Dn40960 &
wait

mv ${data_silo}/S1D_*.txt ${data_ascii}

cd analysis
sed -i -e "s/^\/\/\#define EXCLUDE_HD_MODULE/\#define EXCLUDE_HD_MODULE/" ../../../../source/defines/functionality_flags.hmake -j8

./Radii ${data_silo} S1D_n00080.0 . S1Dn00080
./Radii ${data_silo} S1D_n00160.0 . S1Dn00160
./Radii ${data_silo} S1D_n00320.0 . S1Dn00320
./Radii ${data_silo} S1D_n00640.0 . S1Dn00640
./Radii ${data_silo} S1D_n01280.0 . S1Dn01280
./Radii ${data_silo} S1D_n02560.0 . S1Dn02560
./Radii ${data_silo} S1D_n05120.0 . S1Dn05120
./Radii ${data_silo} S1D_n10240.0 . S1Dn10240
./Radii ${data_silo} S1D_n20480.0 . S1Dn20480
./Radii ${data_silo} S1D_n40960.0 . S1Dn40960

bash plot.sh

cd -
bash plot_density.sh

exit

