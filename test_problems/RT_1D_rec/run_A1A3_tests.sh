#!/bin/bash

test_dir=${1}/RT_1D_rec
code_dir=$2
DATA=$3/RT_1D_rec
mkdir $DATA

cd ${code_dir}
echo "MAKE IN" $code_dir
# compile the code
./compile_code.sh
if [ ! -f ../bin/main_serial ] || [ ! -f ../bin/icgen_serial ]
then
  echo "Cannot compile code"
  exit
else
  echo "MAKE SUCEEDED"
fi
cd $test_dir
rm icgen_serial main_serial
bin_dir=${code_dir}/../bin
cp ${bin_dir}/icgen_serial $test_dir/icgen_serial
cp ${bin_dir}/main_serial $test_dir/main_serial


A1BASE01=SS1D_nh1_dT01_A1
A1BASE02=SS1D_nh1_dT03_A1
A1BASE03=SS1D_nh1_dT10_A1
A1BASE04=SS1D_nh1_dT30_A1
A1BASE11=SS1D_nh2_dT01_A1
A1BASE12=SS1D_nh2_dT03_A1
A1BASE13=SS1D_nh2_dT10_A1
A1BASE14=SS1D_nh2_dT30_A1
A1BASE21=SS1D_nh3_dT01_A1
A1BASE22=SS1D_nh3_dT03_A1
A1BASE23=SS1D_nh3_dT10_A1
A1BASE24=SS1D_nh3_dT30_A1

### TEMP TEMP
################
# make plots   #
################
### TEMP TEMP


TSTEP=( [0]=dt00 [1]=dt01 [2]=dt02 [3]=dt03 [4]=dt04 [5]=dt05 [6]=dt06 [7]=dt07 [8]=dt08 [9]=dt09 [10]=dt10 [11]=dt11 [12]=dt12 )


rm ${DATA}/*SS1D*A1*dt05*

# Generate initial conditions
./icgen_serial pf_${A1BASE01}.txt silo redirect=msg_temp
./icgen_serial pf_${A1BASE02}.txt silo redirect=msg_temp
./icgen_serial pf_${A1BASE03}.txt silo redirect=msg_temp
./icgen_serial pf_${A1BASE04}.txt silo redirect=msg_temp
#
./icgen_serial pf_${A1BASE11}.txt silo redirect=msg_temp
./icgen_serial pf_${A1BASE12}.txt silo redirect=msg_temp
./icgen_serial pf_${A1BASE13}.txt silo redirect=msg_temp
./icgen_serial pf_${A1BASE14}.txt silo redirect=msg_temp
#
./icgen_serial pf_${A1BASE21}.txt silo redirect=msg_temp
./icgen_serial pf_${A1BASE22}.txt silo redirect=msg_temp
./icgen_serial pf_${A1BASE23}.txt silo redirect=msg_temp
./icgen_serial pf_${A1BASE24}.txt silo redirect=msg_temp
#
rm msg_temp*

##############################################
# run models with A1-dt05. #
##############################################
EXE=./main_serial
TSTEP=( [0]=dt00 [1]=dt01 [2]=dt02 [3]=dt03 [4]=dt04 [5]=dt05 [6]=dt06 [7]=dt07 [8]=dt08 [9]=dt09 [10]=dt10 [11]=dt11 [12]=dt12 )


for ii in 5
do
  OPF=1
  if [ $ii -eq 0 ];  then OPF=1;  fi
  if [ $ii -eq 1 ];  then OPF=1;  fi
  if [ $ii -eq 2 ];  then OPF=1;  fi
  if [ $ii -eq 3 ];  then OPF=1;  fi
  if [ $ii -eq 4 ];  then OPF=5;  fi
  if [ $ii -eq 5 ];  then OPF=10; fi
  if [ $ii -eq 6 ];  then OPF=20; fi
  if [ $ii -eq 7 ];  then OPF=40; fi
  if [ $ii -eq 8 ];  then OPF=80; fi
  if [ $ii -eq 9 ];  then OPF=10; fi
  if [ $ii -eq 10 ]; then OPF=20; fi
  if [ $ii -eq 11 ]; then OPF=40; fi
  if [ $ii -eq 12 ]; then OPF=80; fi

  TS=''
  if [ $ii -ge 5 ]; then TS='limit_timestep=5'; fi

  #
  ${EXE} IC_${A1BASE01}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE01}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE01}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE} IC_${A1BASE11}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE11}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE11}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE} IC_${A1BASE21}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE21}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE21}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  wait
  #
  ${EXE} IC_${A1BASE02}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE02}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE02}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE} IC_${A1BASE12}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE12}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE12}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE} IC_${A1BASE22}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE22}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE22}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  wait
  #
  ${EXE} IC_${A1BASE03}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE03}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE03}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE} IC_${A1BASE13}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE13}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE13}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE} IC_${A1BASE23}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE23}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE23}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  wait
  #
  ${EXE} IC_${A1BASE04}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE04}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE04}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE} IC_${A1BASE14}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE14}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE14}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE} IC_${A1BASE24}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE24}_${TSTEP[ii]} redirect=${DATA}/msg_${A1BASE24}_${TSTEP[ii]}_ $TS optype=5 opfreq=$OPF checkpt_freq=100000 &
  wait
done
#exit

############################################
# calculate IF radius as function of time. #
############################################
#
for ii in 5
do
  OPF=1
  if [ $ii -eq 0 ]; then OPF=1;  fi
  if [ $ii -eq 1 ]; then OPF=1;  fi
  if [ $ii -eq 2 ]; then OPF=1;  fi
  if [ $ii -eq 3 ]; then OPF=3;  fi
  if [ $ii -eq 4 ]; then OPF=30; fi
  if [ $ii -eq 5 ]; then OPF=10;  fi
  if [ $ii -eq 6 ]; then OPF=20; fi
  if [ $ii -eq 7 ]; then OPF=40; fi
  if [ $ii -eq 8 ]; then OPF=80; fi
  if [ $ii -eq 9 ]; then OPF=10; fi
  if [ $ii -eq 10 ]; then OPF=20; fi
  if [ $ii -eq 11 ]; then OPF=40; fi
  if [ $ii -eq 12 ]; then OPF=80; fi
  #
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE01}_${TSTEP[ii]} ${DATA}/${A1BASE01}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE02}_${TSTEP[ii]} ${DATA}/${A1BASE02}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE03}_${TSTEP[ii]} ${DATA}/${A1BASE03}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE04}_${TSTEP[ii]} ${DATA}/${A1BASE04}_${TSTEP[ii]}  0 $OPF 5 silo
  #
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE11}_${TSTEP[ii]} ${DATA}/${A1BASE11}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE12}_${TSTEP[ii]} ${DATA}/${A1BASE12}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE13}_${TSTEP[ii]} ${DATA}/${A1BASE13}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE14}_${TSTEP[ii]} ${DATA}/${A1BASE14}_${TSTEP[ii]}  0 $OPF 5 silo
  #
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE21}_${TSTEP[ii]} ${DATA}/${A1BASE21}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE22}_${TSTEP[ii]} ${DATA}/${A1BASE22}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE23}_${TSTEP[ii]} ${DATA}/${A1BASE23}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A1BASE24}_${TSTEP[ii]} ${DATA}/${A1BASE24}_${TSTEP[ii]}  0 $OPF 5 silo
  #
done


##################################################################
############## -------------------------------- ##################
##############  NOW DO THE SAME AGAIN FOR A3  ##################
############## -------------------------------- ##################
##################################################################
A3BASE01=SS1D_nh1_dT01_A3
A3BASE02=SS1D_nh1_dT03_A3
A3BASE03=SS1D_nh1_dT10_A3
A3BASE04=SS1D_nh1_dT30_A3
A3BASE11=SS1D_nh2_dT01_A3
A3BASE12=SS1D_nh2_dT03_A3
A3BASE13=SS1D_nh2_dT10_A3
A3BASE14=SS1D_nh2_dT30_A3
A3BASE21=SS1D_nh3_dT01_A3
A3BASE22=SS1D_nh3_dT03_A3
A3BASE23=SS1D_nh3_dT10_A3
A3BASE24=SS1D_nh3_dT30_A3
TSTEP=( [0]=dt00 [1]=dt01 [2]=dt02 [3]=dt03 [4]=dt04 [5]=dt05 [6]=dt06 [7]=dt07 [8]=dt08 [9]=dt09 [10]=dt10 [11]=dt11 [12]=dt12 )

echo "rm ${DATA}/*A3*dt02*"
rm ${DATA}/*SS1D*A3*dt02*

# Generate initial conditions
./icgen_serial pf_${A3BASE01}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE02}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE03}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE04}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE11}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE12}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE13}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE14}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE21}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE22}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE23}.txt silo redirect=msg_temp
./icgen_serial pf_${A3BASE24}.txt silo redirect=msg_temp
#
rm msg_temp*

##############################################
# run models with A3-dt02. #
##############################################
EXE=./main_serial
TSTEP=( [0]=dt00 [1]=dt01 [2]=dt02 [3]=dt03 [4]=dt04 [5]=dt05 [6]=dt06 [7]=dt07 [8]=dt08 [9]=dt09 [10]=dt10 [11]=dt11 [12]=dt12 )


for ii in 2
do
  OPF=1
  if [ $ii -eq 0 ];  then OPF=4;  fi
  if [ $ii -eq 1 ];  then OPF=8;  fi
  if [ $ii -eq 2 ];  then OPF=16; fi
  if [ $ii -eq 3 ];  then OPF=32; fi
  if [ $ii -eq 4 ];  then OPF=64; fi
  if [ $ii -eq 5 ];  then OPF=20; fi
  if [ $ii -eq 6 ];  then OPF=40; fi
  if [ $ii -eq 7 ];  then OPF=80; fi
  if [ $ii -eq 8 ];  then OPF=160; fi
  if [ $ii -eq 9 ];  then OPF=20; fi
  if [ $ii -eq 10 ]; then OPF=40; fi
  if [ $ii -eq 11 ]; then OPF=80; fi
  if [ $ii -eq 12 ]; then OPF=160; fi
  #
  ${EXE[ii]} IC_${A3BASE01}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE01}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE01}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE[ii]} IC_${A3BASE11}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE11}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE11}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE[ii]} IC_${A3BASE21}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE21}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE21}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  wait
  #
  ${EXE[ii]} IC_${A3BASE02}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE02}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE02}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE[ii]} IC_${A3BASE12}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE12}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE12}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE[ii]} IC_${A3BASE22}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE22}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE22}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  wait
  #
  ${EXE[ii]} IC_${A3BASE03}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE03}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE03}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE[ii]} IC_${A3BASE13}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE13}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE13}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE[ii]} IC_${A3BASE23}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE23}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE23}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  wait
  #
  ${EXE[ii]} IC_${A3BASE04}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE04}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE04}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE[ii]} IC_${A3BASE14}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE14}_${TSTEP[ii]} redirect=${DATA}/msg_${A3BASE14}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  ${EXE[ii]} IC_${A3A3BASE24}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3A3BASE24}_${TSTEP[ii]} redirect=${DATA}/msg_${A3A3BASE24}_${TSTEP[ii]}_ optype=5 opfreq=$OPF checkpt_freq=100000 &
  wait
done
#exit

############################################
# calculate IF radius as function of time. #
############################################
# 
for ii in 2
do
  OPF=1
  if [ $ii -eq 0 ];  then OPF=4;  fi
  if [ $ii -eq 1 ];  then OPF=8;  fi
  if [ $ii -eq 2 ];  then OPF=16; fi
  if [ $ii -eq 3 ];  then OPF=32; fi
  if [ $ii -eq 4 ];  then OPF=64; fi
  if [ $ii -eq 5 ];  then OPF=20; fi
  if [ $ii -eq 6 ];  then OPF=40; fi
  if [ $ii -eq 7 ];  then OPF=80; fi
  if [ $ii -eq 8 ];  then OPF=160; fi
  if [ $ii -eq 9 ];  then OPF=20; fi
  if [ $ii -eq 10 ]; then OPF=40; fi
  if [ $ii -eq 11 ]; then OPF=80; fi
  if [ $ii -eq 12 ]; then OPF=160; fi
  #
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE01}_${TSTEP[ii]} ${DATA}/${A3BASE01}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE02}_${TSTEP[ii]} ${DATA}/${A3BASE02}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE03}_${TSTEP[ii]} ${DATA}/${A3BASE03}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE04}_${TSTEP[ii]} ${DATA}/${A3BASE04}_${TSTEP[ii]}  0 $OPF 5 silo
  #
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE11}_${TSTEP[ii]} ${DATA}/${A3BASE11}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE12}_${TSTEP[ii]} ${DATA}/${A3BASE12}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE13}_${TSTEP[ii]} ${DATA}/${A3BASE13}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE14}_${TSTEP[ii]} ${DATA}/${A3BASE14}_${TSTEP[ii]}  0 $OPF 5 silo
  #
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE21}_${TSTEP[ii]} ${DATA}/${A3BASE21}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE22}_${TSTEP[ii]} ${DATA}/${A3BASE22}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE23}_${TSTEP[ii]} ${DATA}/${A3BASE23}_${TSTEP[ii]}  0 $OPF 5 silo
  ../RT_1D_norec/plot_radius ${DATA}/rad_${A3BASE24}_${TSTEP[ii]} ${DATA}/${A3BASE24}_${TSTEP[ii]}  0 $OPF 5 silo
  #
done

################
# make plots   #
################
echo "PLOTTING PLOTTING PLOTTING PLOTTING V1"
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE01} ${DATA}/rad_${A1BASE01} 6.0e19 3.861e11 3840
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE02} ${DATA}/rad_${A1BASE02} 6.0e19 3.861e11 1280
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE03} ${DATA}/rad_${A1BASE03} 6.0e19 3.861e11 384
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE04} ${DATA}/rad_${A1BASE04} 6.0e19 3.861e11 128
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE11} ${DATA}/rad_${A1BASE11} 6.0e18 3.861e10 3840
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE12} ${DATA}/rad_${A1BASE12} 6.0e18 3.861e10 1280
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE13} ${DATA}/rad_${A1BASE13} 6.0e18 3.861e10 384
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE14} ${DATA}/rad_${A1BASE14} 6.0e18 3.861e10 128
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE21} ${DATA}/rad_${A1BASE21} 6.0e17 3.861e09 3840
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE22} ${DATA}/rad_${A1BASE22} 6.0e17 3.861e09 1280
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE23} ${DATA}/rad_${A1BASE23} 6.0e17 3.861e09 384
./make_A1A3_fig.sh ${DATA}/rad_${A3BASE24} ${DATA}/rad_${A1BASE24} 6.0e17 3.861e09 128

mv *.jpeg *.eps ../

exit











