#!/bin/bash


#OPDIR=/mnt/data/jm/nested_pion/Wind2Dhp
OPDIR=/vol/aibn129/data1/jmackey/scratch/Wind2Dhp
mkdir -p ${OPDIR}

../../icgen_nest_serial params_Wind2Dhp_n0128_l2.txt silo
../../icgen_nest_serial params_Wind2Dhp_n0128_l3.txt silo
../../icgen_nest_serial params_Wind2Dhp_n0128_l4.txt silo
../../icgen_nest_serial params_Wind2Dhp_n0128_l5.txt silo

../../pion_nest_serial fastwind_hp_n0128_l2_level00.00000000.silo outfile=${OPDIR}/OA2_hp_n0128_l2 ooa=2 opfreq=128  redirect=log_OA2_hp_n0128_l2_ &
../../pion_nest_serial fastwind_hp_n0128_l3_level00.00000000.silo outfile=${OPDIR}/OA2_hp_n0128_l3 ooa=2 opfreq=256  redirect=log_OA2_hp_n0128_l3_ &
../../pion_nest_serial fastwind_hp_n0128_l4_level00.00000000.silo outfile=${OPDIR}/OA2_hp_n0128_l4 ooa=2 opfreq=256  redirect=log_OA2_hp_n0128_l4_ &
../../pion_nest_serial fastwind_hp_n0128_l5_level00.00000000.silo outfile=${OPDIR}/OA2_hp_n0128_l5 ooa=2 opfreq=256 redirect=log_OA2_hp_n0128_l5_ &

#../../pion_nest_serial fastwind_hp_n0128_l2_level00.00000000.silo outfile=${OPDIR}/OA1_hp_n0128_l2 ooa=1 opfreq=128  redirect=log_OA1_hp_n0128_l2_ &
#../../pion_nest_serial fastwind_hp_n0128_l3_level00.00000000.silo outfile=${OPDIR}/OA1_hp_n0128_l3 ooa=1 opfreq=256  redirect=log_OA1_hp_n0128_l3_ &
#../../pion_nest_serial fastwind_hp_n0128_l4_level00.00000000.silo outfile=${OPDIR}/OA1_hp_n0128_l4 ooa=1 opfreq=512  redirect=log_OA1_hp_n0128_l4_ &
#../../pion_nest_serial fastwind_hp_n0128_l5_level00.00000000.silo outfile=${OPDIR}/OA1_hp_n0128_l5 ooa=1 opfreq=1024 redirect=log_OA1_hp_n0128_l5_ &

wait

../../icgen_nest_serial params_Wind2Dhp_n0256_l2.txt silo
../../icgen_nest_serial params_Wind2Dhp_n0256_l3.txt silo
../../icgen_nest_serial params_Wind2Dhp_n0256_l4.txt silo
../../icgen_nest_serial params_Wind2Dhp_n0256_l5.txt silo

../../pion_nest_serial fastwind_hp_n0256_l2_level00.00000000.silo outfile=${OPDIR}/OA2_hp_n0256_l2 ooa=2 opfreq=256  redirect=log_OA2_hp_n0256_l2_ &
../../pion_nest_serial fastwind_hp_n0256_l3_level00.00000000.silo outfile=${OPDIR}/OA2_hp_n0256_l3 ooa=2 opfreq=512  redirect=log_OA2_hp_n0256_l3_ &
../../pion_nest_serial fastwind_hp_n0256_l4_level00.00000000.silo outfile=${OPDIR}/OA2_hp_n0256_l4 ooa=2 opfreq=512  redirect=log_OA2_hp_n0256_l4_ &
../../pion_nest_serial fastwind_hp_n0256_l5_level00.00000000.silo outfile=${OPDIR}/OA2_hp_n0256_l5 ooa=2 opfreq=512 redirect=log_OA2_hp_n0256_l5_ &

#../../pion_nest_serial fastwind_hp_n0256_l2_level00.00000000.silo outfile=${OPDIR}/OA1_hp_n0256_l2 ooa=1 opfreq=128  redirect=log_OA1_hp_n0256_l2_ &
#../../pion_nest_serial fastwind_hp_n0256_l3_level00.00000000.silo outfile=${OPDIR}/OA1_hp_n0256_l3 ooa=1 opfreq=256  redirect=log_OA1_hp_n0256_l3_ &
#../../pion_nest_serial fastwind_hp_n0256_l4_level00.00000000.silo outfile=${OPDIR}/OA1_hp_n0256_l4 ooa=1 opfreq=512  redirect=log_OA1_hp_n0256_l4_ &
#../../pion_nest_serial fastwind_hp_n0256_l5_level00.00000000.silo outfile=${OPDIR}/OA1_hp_n0256_l5 ooa=1 opfreq=1024 redirect=log_OA1_hp_n0256_l5_ &

wait

exit


