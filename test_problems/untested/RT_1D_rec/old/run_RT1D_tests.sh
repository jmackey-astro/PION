#!/bin/bash
#
# 2011.04.21 JM: New file, copied from compare_RT_methods.sh in test_RT_nodyn/
# 2011.06.22 JM: Expanded to test the 2nd order convergence of the RT for different MP timestepping.
# 2011.07.09 JM: Adapted file to work on 1D grid.

# current dir
test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/RT_1D_rec
# where exes are located
exe_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/RT_1D_norec
# where compilation happens (not needed here)
code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
# where data files are stored.
data_dir=./data

cd $test_dir
rm ${data_dir}/${fbase0}*

ICGEN=${exe_dir}/icgen_serial
EXE=${exe_dir}/main_serial

###############################################################################
## 1D runs for stromgren spheres with no dynamics, testing the raytracer
## and the energetics (heating and cooling due to H, He, metals).
##  2011.07.09
###############################################################################

fbase0=rtt_Sph1D_ERG_n128_nh2


##
## Generate initial conditions files (uniform neutral medium).
##
$ICGEN ${test_dir}/pf_${fbase0}.txt silo redirect=msg_temp_

##
## Run test problems with recombinations, with different xdot timestep criteria
## First in 1st order accuracy
##
OOA=O1
ACC=${OOA}_DX10000
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

ACC=${OOA}_DX05000
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

ACC=${OOA}_DX02500
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

ACC=${OOA}_DX01250
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

ACC=${OOA}_DX00625
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

# Second order accuracy:
OOA=O2
ACC=${OOA}_DX10000
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

ACC=${OOA}_DX05000
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

ACC=${OOA}_DX02500
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

ACC=${OOA}_DX01250
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

ACC=${OOA}_DX00625
EXE=${exe_dir}/main_serial_${ACC}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 optype=6 checkpt_freq=100000 \
  outfile=${data_dir}/${fbase0}_${ACC} redirect=${data_dir}/msg_${fbase0}_${ACC}_ &

wait

########################################
######### DONE WITH SIMULATIONS ########
########################################
###----------------------------------------------------------------
########################################
echo moving on to analysis

OOA=O1
ACC=${OOA}_DX10000
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo
ACC=${OOA}_DX05000
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo
ACC=${OOA}_DX02500
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo
ACC=${OOA}_DX01250
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo
ACC=${OOA}_DX00625
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo

OOA=O2
ACC=${OOA}_DX10000
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo
ACC=${OOA}_DX05000
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo
ACC=${OOA}_DX02500
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo
ACC=${OOA}_DX01250
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo
ACC=${OOA}_DX00625
${exe_dir}/plot_radius ${data_dir}/${fbase0}_${ACC}  ${data_dir}/${fbase0}_${ACC}  0 5 5 silo


#
#
echo "FINISHED 2D W/ ANALYSIS OF *NEW* 2D SIMS W/ AND W/O RECOMBINATIONS"
#
echo moving on to generate figs comparing different timestepping criteria.

OOA=O2
ACC=${OOA}_DX10000
file1=${data_dir}/${fbase0}_${ACC}
ACC=${OOA}_DX05000
file2=${data_dir}/${fbase0}_${ACC}
ACC=${OOA}_DX02500
file3=${data_dir}/${fbase0}_${ACC}
ACC=${OOA}_DX01250
file4=${data_dir}/${fbase0}_${ACC}
ACC=${OOA}_DX00625
file5=${data_dir}/${fbase0}_${ACC}

FNAME=${test_dir}/COMP1D_radius_error_nh2_n128_${OOA}
cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Ratio R_{a}/R_{if}(t)"
set yrange [0.97:1.015]
set log x
set xrange [0.001:30]
set title ""
set output "${FNAME}.eps"
set label "1D, n_H=100 cm^{-3}, 2nd Order" at 0.25,1.0075
set label "Cell {/Symbol d}{/Symbol t}=30 for N=128 runs" at 0.25,1.0125
g(x)=(1.0-exp(-x))**(1.0/3.0)
RH=6.0e18
RI=4.5172e18
tH=3.861e10
tI=3.861e10
f(x)=1.0
h(x)=1.0+4.6875e16/(RI*g(x))
j(x)=1.0-4.6875e16/(RI*g(x))
plot f(x) w l lw 2 lt -1 title 'Correct Solution', h(x) w l lw 2 lt 0 title "One cell error", j(x) w l lw 2 lt 0 notitle, \
  '${file4}.txt' u (\$1/tI):((\$5/RI)/g(\$1/tI)) w l lw 2 lt 2  title 'N=128, {/Symbol d}t= 0.125/xdot', \
  '${file3}.txt' u (\$1/tI):((\$5/RI)/g(\$1/tI)) w l lw 2 lt 3  title 'N=128, {/Symbol d}t= 0.250/xdot', \
  '${file2}.txt' u (\$1/tI):((\$5/RI)/g(\$1/tI)) w l lw 2 lt 1  title 'N=128, {/Symbol d}t= 0.500/xdot'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg


# FIRST ORDER ACCURACY (LARGER RANGE IN Y-AXIS)
OOA=O1
ACC=${OOA}_DX10000
file1=${data_dir}/${fbase0}_${ACC}
ACC=${OOA}_DX05000
file2=${data_dir}/${fbase0}_${ACC}
ACC=${OOA}_DX02500
file3=${data_dir}/${fbase0}_${ACC}
ACC=${OOA}_DX01250
file4=${data_dir}/${fbase0}_${ACC}
ACC=${OOA}_DX00625
file5=${data_dir}/${fbase0}_${ACC}

FNAME=${test_dir}/COMP1D_radius_error_nh2_n128_${OOA}
cat << EOF  > gnu.plt
set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Ratio R_{a}/R_{if}(t)"
set yrange [0.89:1.02]
set log x
set xrange [0.001:30]
set title ""
set output "${FNAME}.eps"
set label "1D, n_H=100 cm^{-3}, 1st Order" at 0.0015,1.0075
set label "Cell {/Symbol d}{/Symbol t}=30 for N=128 runs" at 0.0015,1.015
g(x)=(1.0-exp(-x))**(1.0/3.0)
RH=6.0e18
RI=4.5172e18
tH=3.861e10
tI=3.861e10
f(x)=1.0
h(x)=1.0+4.6875e16/(RI*g(x))
j(x)=1.0-4.6875e16/(RI*g(x))
plot f(x) w l lw 2 lt -1 title 'Correct Solution', h(x) w l lw 2 lt 0 title "One cell error", j(x) w l lw 2 lt 0 notitle, \
  '${file4}.txt' u (\$1/tI):((\$5/RI)/g(\$1/tI)) w l lw 2 lt 2  title 'N=128, {/Symbol d}t= 0.125/xdot', \
  '${file3}.txt' u (\$1/tI):((\$5/RI)/g(\$1/tI)) w l lw 2 lt 3  title 'N=128, {/Symbol d}t= 0.250/xdot', \
  '${file2}.txt' u (\$1/tI):((\$5/RI)/g(\$1/tI)) w l lw 2 lt 1  title 'N=128, {/Symbol d}t= 0.500/xdot'
#pause -1
EOF
gnuplot gnu.plt
convert -density 300 -quality 100 ${FNAME}.eps ${FNAME}.jpeg


rm ${data_dir}/${fbase0}*

exit

###----------------------------------------------------------------


