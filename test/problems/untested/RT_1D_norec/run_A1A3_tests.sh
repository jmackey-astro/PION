#!/bin/bash

# 2011.11.25 JM: Generate Fig1 for the RT-algorithms paper.

test_dir=${1}/RT_1D_norec
code_dir=$2
DATA=$3/RT_1D_norec
mkdir $DATA

A1BASE00=rtt_NoRec_cart1D_dT10_v0010_n1e1_ALG1
A1BASE01=rtt_NoRec_cart1D_dT10_v0010_n1e2_ALG1
A1BASE02=rtt_NoRec_cart1D_dT10_v0010_n1e3_ALG1
A1BASE03=rtt_NoRec_cart1D_dT10_v0030_n1e1_ALG1
A1BASE04=rtt_NoRec_cart1D_dT10_v0030_n1e2_ALG1
A1BASE05=rtt_NoRec_cart1D_dT10_v0030_n1e3_ALG1
A1BASE06=rtt_NoRec_cart1D_dT10_v0100_n1e1_ALG1
A1BASE07=rtt_NoRec_cart1D_dT10_v0100_n1e2_ALG1
A1BASE08=rtt_NoRec_cart1D_dT10_v0100_n1e3_ALG1
A1BASE09=rtt_NoRec_cart1D_dT10_v0300_n1e1_ALG1
A1BASE10=rtt_NoRec_cart1D_dT10_v0300_n1e2_ALG1
A1BASE11=rtt_NoRec_cart1D_dT10_v0300_n1e3_ALG1
A1BASE12=rtt_NoRec_cart1D_dT10_v1000_n1e1_ALG1
A1BASE13=rtt_NoRec_cart1D_dT10_v1000_n1e2_ALG1
A1BASE14=rtt_NoRec_cart1D_dT10_v1000_n1e3_ALG1

A3BASE00=rtt_NoRec_cart1D_dT10_v0010_n1e1_ALG3
A3BASE01=rtt_NoRec_cart1D_dT10_v0010_n1e2_ALG3
A3BASE02=rtt_NoRec_cart1D_dT10_v0010_n1e3_ALG3
A3BASE03=rtt_NoRec_cart1D_dT10_v0030_n1e1_ALG3
A3BASE04=rtt_NoRec_cart1D_dT10_v0030_n1e2_ALG3
A3BASE05=rtt_NoRec_cart1D_dT10_v0030_n1e3_ALG3
A3BASE06=rtt_NoRec_cart1D_dT10_v0100_n1e1_ALG3
A3BASE07=rtt_NoRec_cart1D_dT10_v0100_n1e2_ALG3
A3BASE08=rtt_NoRec_cart1D_dT10_v0100_n1e3_ALG3
A3BASE09=rtt_NoRec_cart1D_dT10_v0300_n1e1_ALG3
A3BASE10=rtt_NoRec_cart1D_dT10_v0300_n1e2_ALG3
A3BASE11=rtt_NoRec_cart1D_dT10_v0300_n1e3_ALG3
A3BASE12=rtt_NoRec_cart1D_dT10_v1000_n1e1_ALG3
A3BASE13=rtt_NoRec_cart1D_dT10_v1000_n1e2_ALG3
A3BASE14=rtt_NoRec_cart1D_dT10_v1000_n1e3_ALG3

####### TEMP ###########
####### TEMP ###########


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

cd ../test_RT_nodyn/
make -f Makefile.plotradius
cp plot_radius $test_dir
cd $test_dir


###########################
### ALG1 SET FILE NAMES ###
###########################
A1BASE=( [0]=${A1BASE00} [1]=${A1BASE01} [2]=${A1BASE02} [3]=${A1BASE03} [4]=${A1BASE04} [5]=${A1BASE05} [6]=${A1BASE06} [7]=${A1BASE07} [8]=${A1BASE08} [9]=${A1BASE09} [10]=${A1BASE10} [11]=${A1BASE11} [12]=${A1BASE12} [13]=${A1BASE13} [14]=${A1BASE14} )
TSTEP=( [0]=dt00 [1]=dt01 [2]=dt02 [3]=dt03 [4]=dt04 [5]=dt05 [6]=dt06 [7]=dt07 [8]=dt08 [9]=dt09 [10]=dt10 [11]=dt11 [12]=dt12 )
EXE=./main_serial

## delete any pre-existing files
for ff in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
  rm ${DATA}/*${A1BASE[ff]}*
done

## GENERATE INITIAL CONDITIONS
for ff in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
  ./icgen_serial pf_${A1BASE[ff]}.txt silo redirect=msg_temp_
done
rm msg_temp_*

## LOOP OVER dt05
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
  ## LOOP OVER ALL SIMULATIONS
  for ff in 0 3 6 9 12
  do
    let f0=$ff
    let f1=$ff+1
    let f2=$ff+2
    ${EXE} IC_${A1BASE[f0]}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE[f0]}_${TSTEP[ii]} \
      redirect=${DATA}/msg_${A1BASE[f0]}_${TSTEP[ii]}_ \
      optype=5 opfreq=$OPF $TS checkpt_freq=100000 &
    ${EXE} IC_${A1BASE[f1]}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE[f1]}_${TSTEP[ii]} \
      redirect=${DATA}/msg_${A1BASE[f1]}_${TSTEP[ii]}_ \
      optype=5 opfreq=$OPF $TS checkpt_freq=100000 &
    ${EXE} IC_${A1BASE[f2]}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A1BASE[f2]}_${TSTEP[ii]} \
      redirect=${DATA}/msg_${A1BASE[f2]}_${TSTEP[ii]}_ \
      optype=5 opfreq=$OPF $TS checkpt_freq=100000 &
    wait
    ./plot_radius ${DATA}/results_${A1BASE[f0]}_${TSTEP[ii]} ${DATA}/${A1BASE[f0]}_${TSTEP[ii]} 0 $OPF 5 silo
    ./plot_radius ${DATA}/results_${A1BASE[f1]}_${TSTEP[ii]} ${DATA}/${A1BASE[f1]}_${TSTEP[ii]} 0 $OPF 5 silo
    ./plot_radius ${DATA}/results_${A1BASE[f2]}_${TSTEP[ii]} ${DATA}/${A1BASE[f2]}_${TSTEP[ii]} 0 $OPF 5 silo
  done
done


#####################
### ALG3 ###
#####################
A3BASE=( [0]=${A3BASE00} [1]=${A3BASE01} [2]=${A3BASE02} [3]=${A3BASE03} [4]=${A3BASE04} [5]=${A3BASE05} [6]=${A3BASE06} [7]=${A3BASE07} [8]=${A3BASE08} [9]=${A3BASE09} [10]=${A3BASE10} [11]=${A3BASE11} [12]=${A3BASE12} [13]=${A3BASE13} [14]=${A3BASE14} )
EXE=./main_serial
TSTEP=( [0]=dt00 [1]=dt01 [2]=dt02 [3]=dt03 [4]=dt04 [5]=dt05 [6]=dt06 [7]=dt07 [8]=dt08 [9]=dt09 [10]=dt10 [11]=dt11 [12]=dt12 )

for ff in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
  rm ${DATA}/*${A3BASE[ff]}*
done

for ff in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
  ./icgen_serial pf_${A3BASE[ff]}.txt silo redirect=msg_temp_
done
rm msg_temp_*

## RUN ALG3 SIMULATIONS for dt02
for ii in 2
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
  
  ## LOOP OVER ALL SIMULATIONS (LEAVING OUT 12-14!!)
  for ff in 0 3 6 9 12
  do
    let f0=$ff
    let f1=$ff+1
    let f2=$ff+2
    ${EXE} IC_${A3BASE[f0]}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE[f0]}_${TSTEP[ii]} \
      redirect=${DATA}/msg_${A3BASE[f0]}_${TSTEP[ii]}_ \
      optype=5 opfreq=$OPF checkpt_freq=100000 &
    ${EXE} IC_${A3BASE[f1]}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE[f1]}_${TSTEP[ii]} \
      redirect=${DATA}/msg_${A3BASE[f1]}_${TSTEP[ii]}_ \
      optype=5 opfreq=$OPF checkpt_freq=100000 &
    ${EXE} IC_${A3BASE[f2]}.silo 5 1 cfl=10000.0 outfile=${DATA}/${A3BASE[f2]}_${TSTEP[ii]} \
      redirect=${DATA}/msg_${A3BASE[f2]}_${TSTEP[ii]}_ \
      optype=5 opfreq=$OPF checkpt_freq=100000 &
    wait
    ./plot_radius ${DATA}/results_${A3BASE[f0]}_${TSTEP[ii]} ${DATA}/${A3BASE[f0]}_${TSTEP[ii]} 0 $OPF 5 silo
    ./plot_radius ${DATA}/results_${A3BASE[f1]}_${TSTEP[ii]} ${DATA}/${A3BASE[f1]}_${TSTEP[ii]} 0 $OPF 5 silo
    ./plot_radius ${DATA}/results_${A3BASE[f2]}_${TSTEP[ii]} ${DATA}/${A3BASE[f2]}_${TSTEP[ii]} 0 $OPF 5 silo
  done
done




################################################################################
# Plot error for each simulation for A1 and A3 
################################################################################
A1BASE=( [0]=${A1BASE00} [1]=${A1BASE01} [2]=${A1BASE02} [3]=${A1BASE03} [4]=${A1BASE04} [5]=${A1BASE05} [6]=${A1BASE06} [7]=${A1BASE07} [8]=${A1BASE08} [9]=${A1BASE09} [10]=${A1BASE10} [11]=${A1BASE11} [12]=${A1BASE12} [13]=${A1BASE13} [14]=${A1BASE14} )

A3BASE=( [0]=${A3BASE00} [1]=${A3BASE01} [2]=${A3BASE02} [3]=${A3BASE03} [4]=${A3BASE04} [5]=${A3BASE05} [6]=${A3BASE06} [7]=${A3BASE07} [8]=${A3BASE08} [9]=${A3BASE09} [10]=${A3BASE10} [11]=${A3BASE11} [12]=${A3BASE12} [13]=${A3BASE13} [14]=${A3BASE14} )

for ff in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 
do
  FNAME=ACC_${A3BASE[ff]}ALG1
  cat << EOF  > gnu.plt
  set terminal postscript enhanced color eps
  set size 0.7071,0.7071
  set xlabel "Time (kyr)" 0.0,0.5
  set ylabel "Fractional Error" 3.0,1.0
  set log y
  set yrange [1.0e-6:1.0]
  unset log x
  #set xrange [0:12]
  #set xtics 0,1,12
  set title ""
  set output "${FNAME}.eps"
  plot '${DATA}/results_${A1BASE[ff]}_dt05.txt' u (\$1/3.16e10):(abs(1.0-\$6/\$7)) w lp lw 2 title "A1", \
       '${DATA}/results_${A3BASE[ff]}_dt02.txt' u (\$1/3.16e10):(abs(1.0-\$6/\$7)) w lp lw 2 title "A3"
  exit
EOF
  gnuplot gnu.plt
  convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg
done

mv *.jpeg *.eps ../

exit

