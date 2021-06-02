#!/bin/bash

###########
# recombs
###########
cp pf_rtt2d_n32_nh1_rec.txt   pf_rtt2d_n32_nh2_rec.txt  
cp pf_rtt2d_n100_nh1_rec.txt  pf_rtt2d_n100_nh2_rec.txt
cp pf_rtt2d_n256_nh1_rec.txt  pf_rtt2d_n256_nh2_rec.txt  
cp pf_rtt3d_n32_nh1_rec.txt   pf_rtt3d_n32_nh2_rec.txt
cp pf_rtt3d_n100_nh1_rec.txt  pf_rtt3d_n100_nh2_rec.txt  

sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt2d_n32_nh1_rec.txt
sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt2d_n100_nh1_rec.txt
sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt2d_n256_nh1_rec.txt
sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt3d_n32_nh1_rec.txt
sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt3d_n100_nh1_rec.txt
sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt2d_n32_nh2_rec.txt
sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt2d_n100_nh2_rec.txt
sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt2d_n256_nh2_rec.txt
sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt3d_n32_nh2_rec.txt
sed -i -e "s/OutputFileType fits/OutputFileType silo/g" pf_rtt3d_n100_nh2_rec.txt

sed -i -e "s/PEC_ambRO 1.67e-23/PEC_ambRO 1.67e-22/g" pf_rtt2d_n32_nh2_rec.txt  
sed -i -e "s/PEC_ambRO 1.67e-23/PEC_ambRO 1.67e-22/g" pf_rtt2d_n100_nh2_rec.txt  
sed -i -e "s/PEC_ambRO 1.67e-23/PEC_ambRO 1.67e-22/g" pf_rtt2d_n256_nh2_rec.txt  
sed -i -e "s/PEC_ambRO 1.67e-23/PEC_ambRO 1.67e-22/g" pf_rtt3d_n32_nh2_rec.txt  
sed -i -e "s/PEC_ambRO 1.67e-23/PEC_ambRO 1.67e-22/g" pf_rtt3d_n100_nh2_rec.txt  

sed -i -e "s/PEC_ambPG 2.505e-12/PEC_ambPG 2.505e-11/g" pf_rtt2d_n32_nh2_rec.txt
sed -i -e "s/PEC_ambPG 2.505e-12/PEC_ambPG 2.505e-11/g" pf_rtt2d_n100_nh2_rec.txt
sed -i -e "s/PEC_ambPG 2.505e-12/PEC_ambPG 2.505e-11/g" pf_rtt2d_n256_nh2_rec.txt
sed -i -e "s/PEC_ambPG 2.505e-12/PEC_ambPG 2.505e-11/g" pf_rtt3d_n32_nh2_rec.txt
sed -i -e "s/PEC_ambPG 2.505e-12/PEC_ambPG 2.505e-11/g" pf_rtt3d_n100_nh2_rec.txt

sed -i -e "s/nh1_rec/nh2_rec/g" pf_rtt2d_n32_nh2_rec.txt
sed -i -e "s/nh1_rec/nh2_rec/g" pf_rtt2d_n100_nh2_rec.txt
sed -i -e "s/nh1_rec/nh2_rec/g" pf_rtt2d_n256_nh2_rec.txt
sed -i -e "s/nh1_rec/nh2_rec/g" pf_rtt3d_n32_nh2_rec.txt
sed -i -e "s/nh1_rec/nh2_rec/g" pf_rtt3d_n100_nh2_rec.txt

sed -i -e "s/RPstrength0 3.0e25/RPstrength0 3.0e27/g" pf_rtt2d_n32_nh2_rec.txt
sed -i -e "s/RPstrength0 3.0e25/RPstrength0 3.0e27/g" pf_rtt2d_n100_nh2_rec.txt
sed -i -e "s/RPstrength0 3.0e25/RPstrength0 3.0e27/g" pf_rtt2d_n256_nh2_rec.txt
sed -i -e "s/RPstrength0 3.0e43/RPstrength0 3.0e45/g" pf_rtt3d_n32_nh2_rec.txt
sed -i -e "s/RPstrength0 3.0e43/RPstrength0 3.0e45/g" pf_rtt3d_n100_nh2_rec.txt

sed -i -e "s/FinishTime 3.861e12/FinishTime 3.861e11/g" pf_rtt2d_n32_nh2_rec.txt
sed -i -e "s/FinishTime 3.861e12/FinishTime 3.861e11/g" pf_rtt2d_n100_nh2_rec.txt
sed -i -e "s/FinishTime 3.861e12/FinishTime 3.861e11/g" pf_rtt2d_n256_nh2_rec.txt
sed -i -e "s/FinishTime 3.861e12/FinishTime 3.861e11/g" pf_rtt3d_n32_nh2_rec.txt
sed -i -e "s/FinishTime 3.861e12/FinishTime 3.861e11/g" pf_rtt3d_n100_nh2_rec.txt

exit


###########
#
# 3D data
#
cp pf_rtt3d_n32_nh1_norec.txt  pf_rtt3d_n32_nh2_norec.txt
cp pf_rtt3d_n100_nh1_norec.txt pf_rtt3d_n100_nh2_norec.txt
cp pf_rtt3d_n32_nh1_norec.txt  pf_rtt3d_n32_nh3_norec.txt
cp pf_rtt3d_n100_nh1_norec.txt pf_rtt3d_n100_nh3_norec.txt

sed -i -e "s/PEC_ambRO 1.67e-23/PEC_ambRO 1.67e-22/g" pf_rtt3d_n32_nh2_norec.txt
sed -i -e "s/PEC_ambRO 1.67e-23/PEC_ambRO 1.67e-22/g" pf_rtt3d_n100_nh2_norec.txt
sed -i -e "s/PEC_ambRO 1.67e-23/PEC_ambRO 1.67e-21/g" pf_rtt3d_n32_nh3_norec.txt
sed -i -e "s/PEC_ambRO 1.67e-23/PEC_ambRO 1.67e-21/g" pf_rtt3d_n100_nh3_norec.txt

sed -i -e "s/PEC_ambPG 2.505e-12/PEC_ambPG 2.505e-11/g" pf_rtt3d_n32_nh2_norec.txt
sed -i -e "s/PEC_ambPG 2.505e-12/PEC_ambPG 2.505e-11/g" pf_rtt3d_n100_nh2_norec.txt
sed -i -e "s/PEC_ambPG 2.505e-12/PEC_ambPG 2.505e-10/g" pf_rtt3d_n32_nh3_norec.txt
sed -i -e "s/PEC_ambPG 2.505e-12/PEC_ambPG 2.505e-10/g" pf_rtt3d_n100_nh3_norec.txt

sed -i -e "s/nh1_norec/nh2_norec/g" pf_rtt3d_n32_nh2_norec.txt
sed -i -e "s/nh1_norec/nh2_norec/g" pf_rtt3d_n100_nh2_norec.txt
sed -i -e "s/nh1_norec/nh3_norec/g" pf_rtt3d_n32_nh3_norec.txt
sed -i -e "s/nh1_norec/nh3_norec/g" pf_rtt3d_n100_nh3_norec.txt

sed -i -e "s/RPstrength0 3.0e43/RPstrength0 3.0e44/g" pf_rtt3d_n32_nh2_norec.txt
sed -i -e "s/RPstrength0 3.0e43/RPstrength0 3.0e44/g" pf_rtt3d_n100_nh2_norec.txt
sed -i -e "s/RPstrength0 3.0e43/RPstrength0 3.0e45/g" pf_rtt3d_n32_nh3_norec.txt
sed -i -e "s/RPstrength0 3.0e43/RPstrength0 3.0e45/g" pf_rtt3d_n100_nh3_norec.txt


exit

#
# 2D data
#
cp pf_rtt2d_n32_nh2_norec.txt  pf_rtt2d_n32_nh3_norec.txt
cp pf_rtt2d_n100_nh2_norec.txt pf_rtt2d_n100_nh3_norec.txt
cp pf_rtt2d_n256_nh2_norec.txt pf_rtt2d_n256_nh3_norec.txt

sed -i -e "s/PEC_ambRO 1.67e-22/PEC_ambRO 1.67e-21/g" pf_rtt2d_n32_nh3_norec.txt
sed -i -e "s/PEC_ambRO 1.67e-22/PEC_ambRO 1.67e-21/g" pf_rtt2d_n100_nh3_norec.txt
sed -i -e "s/PEC_ambRO 1.67e-22/PEC_ambRO 1.67e-21/g" pf_rtt2d_n256_nh3_norec.txt

sed -i -e "s/PEC_ambPG 2.505e-11/PEC_ambPG 2.505e-10/g" pf_rtt2d_n32_nh3_norec.txt
sed -i -e "s/PEC_ambPG 2.505e-11/PEC_ambPG 2.505e-10/g" pf_rtt2d_n100_nh3_norec.txt
sed -i -e "s/PEC_ambPG 2.505e-11/PEC_ambPG 2.505e-10/g" pf_rtt2d_n256_nh3_norec.txt

sed -i -e "s/nh2_norec/nh3_norec/g" pf_rtt2d_n32_nh3_norec.txt
sed -i -e "s/nh2_norec/nh3_norec/g" pf_rtt2d_n100_nh3_norec.txt
sed -i -e "s/nh2_norec/nh3_norec/g" pf_rtt2d_n256_nh3_norec.txt

sed -i -e "s/RPstrength0 3.0e26/RPstrength0 3.0e27/g" pf_rtt2d_n32_nh3_norec.txt
sed -i -e "s/RPstrength0 3.0e26/RPstrength0 3.0e27/g" pf_rtt2d_n100_nh3_norec.txt
sed -i -e "s/RPstrength0 3.0e26/RPstrength0 3.0e27/g" pf_rtt2d_n256_nh3_norec.txt

sed -i -e "s/Tracer0 1.0e-6/Tracer0 1.0e-12/g" pf_rtt2d_n32_nh3_norec.txt
sed -i -e "s/Tracer0 1.0e-6/Tracer0 1.0e-12/g" pf_rtt2d_n100_nh3_norec.txt
sed -i -e "s/Tracer0 1.0e-6/Tracer0 1.0e-12/g" pf_rtt2d_n256_nh3_norec.txt
