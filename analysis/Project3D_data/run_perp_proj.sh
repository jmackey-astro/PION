#!/bin/bash
#
# 2012.05.11 JM: Running on JUROPA to make plots of moving HII regions.
#

DATA_DIR=/vol/trillian2/raid5/jmackey/MovingHIIregion/ZetaOph
#DATA_DIR=/lustre/jwork/hbn23/hbn231/MovingHIIregion/ZetaOph
OP_DIR=/vol/klaipeda2/scratch2/jmackey/MovingHIIregion_projections
mkdir $OP_DIR
CODE_DIR=.

cd $CODE_DIR
#
# projection <data_dir> <input-file> <output-path+file> <ftype=1(FITS)> 
#            <multi-files=1> <LOS> <VERT> <ANGLE-to-LOS> <What2int>
#            [<Nbins> <v_min> <v_max> <smooth?> <how-much>]
#
#
VLOS=2
HALPHA=4
DENS=0
NEUTRALDENS=1
ALL_SCALARS=7

#
# dr=0.30pc lowest resolution models.
#
#FILE1=ZOph3D_Bpll7uG_dr0p3_NearlyAlignedv2_NoNs_DS1p0em6_RCV_0000.00000
#FITSFILE1=ZOph3D_Bpll7uG_dr0p30_XPYP_ZNlos_NTDENS
#time ./projection ${DATA_DIR} $FILE1 \
# ${OP_DIR}/${FITSFILE1} 1 1 ZN YP 00 $NEUTRALDENS &

#FILE2=ZOph3D_Bpll7uG_dr0p3_NearlyAlignedv2_NoNs_DS1p0em6_RCV_0000.00000
#FITSFILE2=ZOph3D_Bpll7uG_dr0p30_XPYP_ZNlos_DENS
#time ./projection ${DATA_DIR} $FILE2 \
# ${OP_DIR}/${FITSFILE2} 1 1 ZN YP 00 $DENS &

FILE3=ZOph3D_Bpll7uG_dr0p3_NearlyAlignedv2_NoNs_DS1p0em6_RCV_0000.00000
FITSFILE3=ZOph3D_Bpll7uG_dr0p30_XPYP_ZNlos_Emission
time ./projection ${DATA_DIR} $FILE3 \
 ${OP_DIR}/${FITSFILE3} 1 1 ZN YP 00 $HALPHA &

wait

#grep "SIMULATION TIME" ${OP_DIR}/${FITSFILE1}_msg_info.txt > ${OP_DIR}/${FITSFILE1}_timesteps.txt 
#grep "SIMULATION TIME" ${OP_DIR}/${FITSFILE2}_msg_info.txt > ${OP_DIR}/${FITSFILE2}_timesteps.txt
grep "SIMULATION TIME" ${OP_DIR}/${FITSFILE3}_msg_info.txt > ${OP_DIR}/${FITSFILE3}_timesteps.txt

exit

#
# dr=0.30pc lowest resolution models.
#
FILE1=ZOph3D_Hydro_dr0p3_NearlyAlignedv2_NoNs_DS1p0em6_HYB_0000.0000
FITSFILE1=ZOph3D_Hydro_dr0p30_XPYP_ZNlos_NTDENS
time ./projection ${DATA_DIR} $FILE1 \
 ${OP_DIR}/${FITSFILE1} 1 1 ZN YP 00 $NEUTRALDENS &

FILE2=ZOph3D_Hydro_dr0p3_NearlyAlignedv2_NoNs_DS1p0em6_HYB_0000.0000
FITSFILE2=ZOph3D_Hydro_dr0p30_XPYP_ZNlos_DENS
time ./projection ${DATA_DIR} $FILE2 \
 ${OP_DIR}/${FITSFILE2} 1 1 ZN YP 00 $DENS &

FILE3=ZOph3D_Hydro_dr0p3_NearlyAlignedv2_NoNs_DS1p0em6_HYB_0000.0000
FITSFILE3=ZOph3D_Hydro_dr0p30_XPYP_ZNlos_Emission
time ./projection ${DATA_DIR} $FILE3 \
 ${OP_DIR}/${FITSFILE3} 1 1 ZN YP 00 $HALPHA &

wait

grep "SIMULATION TIME" ${OP_DIR}/${FITSFILE1}_msg_info.txt > ${OP_DIR}/${FITSFILE1}_timesteps.txt 
grep "SIMULATION TIME" ${OP_DIR}/${FITSFILE2}_msg_info.txt > ${OP_DIR}/${FITSFILE2}_timesteps.txt
grep "SIMULATION TIME" ${OP_DIR}/${FITSFILE3}_msg_info.txt > ${OP_DIR}/${FITSFILE3}_timesteps.txt

wait

exit

#
# dr=0.15pc medium resolution models.
#
FILE1=ZOph3D_Hydro_dr0p15_NearlyAlignedv2_NoNs_DS1p0em6_HYB_0000.0000
FITSFILE1=ZOph3D_Hydro_dr0p15_XPYP_ZNlos_NTDENS
time ./projection ${DATA_DIR} $FILE1 \
 ${OP_DIR}/${FITSFILE1} 1 1 ZN YP 00 $NEUTRALDENS

FILE2=ZOph3D_Hydro_dr0p15_NearlyAlignedv2_NoNs_DS1p0em6_HYB_0000.0000
FITSFILE2=ZOph3D_Hydro_dr0p15_XPYP_ZNlos_DENS
time ./projection ${DATA_DIR} $FILE2 \
 ${OP_DIR}/${FITSFILE2} 1 1 ZN YP 00 $DENS

FILE3=ZOph3D_Hydro_dr0p15_NearlyAlignedv2_NoNs_DS1p0em6_HYB_0000.0000
FITSFILE3=ZOph3D_Hydro_dr0p15_XPYP_ZNlos_Emission
time ./projection ${DATA_DIR} $FILE3 \
 ${OP_DIR}/${FITSFILE3} 1 1 ZN YP 00 $HALPHA

grep "SIMULATION TIME" ${OP_DIR}/${FITSFILE1}_msg_info.txt > ${OP_DIR}/${FITSFILE1}_timesteps.txt 
grep "SIMULATION TIME" ${OP_DIR}/${FITSFILE2}_msg_info.txt > ${OP_DIR}/${FITSFILE2}_timesteps.txt
grep "SIMULATION TIME" ${OP_DIR}/${FITSFILE3}_msg_info.txt > ${OP_DIR}/${FITSFILE3}_timesteps.txt

exit


