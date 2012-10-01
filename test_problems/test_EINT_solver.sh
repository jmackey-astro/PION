#!/bin/bash
#
# 2010.12.31 JM: Run all hydro tests using the internal/total energy
#  hybrid solver.  Copied from run_all_tests.sh
#




cd ../../
BASE_DIR=`pwd`
cd -
code_dir=${BASE_DIR}/trunk/bin_serial
test_dir=${BASE_DIR}/trunk/test_problems
# If the code is on a network drive, code may run faster with data_dir
# set to a local disk.
data_dir=/vol/aibn129/aibn129_1/jmackey/data_etc/code_tests/EINT2

# cmp_dir doesn't do anything yet.  It will contain some results from
# standard problems to compare the current tests to.
cmp_dir=/vol/aibn129/aibn129_1/jmackey/data_etc/code_tests/ref
visit_cmd=/vol/aibn129/aibn129_1/jmackey/extra_libraries/visit_bin/bin/visit

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir
# check that directory exists? how?

#
# We want to set a testing flag, so the code does not calculate
# boundary-point timesteps in the first 10 timesteps, and also sets a
# fixed cross-section for hydrogen ionisation and recombination.
#
sed -i -e "s|//#define RT_TEST_PROBS|#define RT_TEST_PROBS|g" ${BASE_DIR}/trunk/source/global.h
#
# Check we can compile the code
#
cd $code_dir
echo "MAKE IN" $code_dir
#make -f Makefile.serial.code; make -f Makefile.serial.icgenerator
## check that executables exist? how?
#echo "MAKE SUCEEDED"

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Test Problem: Blastwave (1D,2D (slab+axi),3D)
#
echo "*********** BLAST--WAVE TESTS!!! PUT YOUR GOGGLES ON ;-) **********"
##
## 1D spherical blast wave goes here
##
cd ${test_dir}/blastwave_sph1d
tdir=${test_dir}/blastwave_sph1d
ddir=${data_dir}/blastwave_sph1d
# In case it doesn't exist, create the destination directory.
mkdir $ddir

sed -i -e "s|//#define BLAST_WAVE_CHECK|#define BLAST_WAVE_CHECK|g" ${BASE_DIR}/trunk/source/global.h
cd ${code_dir}
echo "MAKE IN" $code_dir
make -f Makefile.serial.code; make -f Makefile.serial.icgenerator
echo "MAKE SUCEEDED"
./icgen ${tdir}/pf_sphBW_n128.txt   silo redirect=tmp_
./icgen ${tdir}/pf_sphBW_n256.txt   silo redirect=tmp_
./icgen ${tdir}/pf_sphBW_n512.txt   silo redirect=tmp_

#
# N=128
./main_serial IC_BW1D_phys_n128.silo 5 1 outfile=${ddir}/BW1D_phys_n128_EINT_FKJav01 \
 redirect=${ddir}/msg_BW1D_phys_n128_EINT_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=3 eqntype=9 &
#
# N=256
./main_serial IC_BW1D_phys_n256.silo 5 1 outfile=${ddir}/BW1D_phys_n256_EINT_FKJav01 \
 redirect=${ddir}/msg_BW1D_phys_n256_EINT_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=3 eqntype=9 &
#
# N=512
./main_serial IC_BW1D_phys_n512.silo 5 1 outfile=${ddir}/BW1D_phys_n512_EINT_FKJav01 \
 redirect=${ddir}/msg_BW1D_phys_n512_EINT_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=3 eqntype=9 &
#
wait
cd $tdir

file1=`ls ${ddir}/BW1D_phys_n128_EINT_FKJav01.0* | tail -n1`
file2=`ls ${ddir}/BW1D_phys_n256_EINT_FKJav01.0* | tail -n1`
file3=`ls ${ddir}/BW1D_phys_n512_EINT_FKJav01.0* | tail -n1`
./plot_sphBW_dens_pres_profile.sh BW1D_Profile50kyr_EINT_FKJ98 \
  $file1 $file2 $file3 \
  "EINT, FKJ98, N=128" "EINT, FKJ98, N=256" "EINT, FKJ98, N=512"

echo "*********** 1D BW TESTS: FIGS GENERATED (SEE blastwave_sph1d/*.jpeg) **********"
eog ${tdir}/*EINT*.jpeg &
sed -i -e "s|#define BLAST_WAVE_CHECK|//#define BLAST_WAVE_CHECK|g" ${BASE_DIR}/trunk/source/global.h
#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Adiabatic 2D,3D blast wave tests here
#---------------------------------------------------------------------
cd ${test_dir}/blastwave_axi2d
tdir=${test_dir}/blastwave_axi2d
ddir=${data_dir}/blastwave_axi2d
# In case it doesn't exist, create the destination directory.
mkdir $ddir
cd ${code_dir}
echo "MAKE IN" $code_dir
make -f Makefile.serial.code; make -f Makefile.serial.icgenerator
echo "MAKE SUCEEDED"
#
./icgen ${tdir}/pf_axi2dBW_HalfPlane_NR064.txt   silo redirect=tmp_
./icgen ${tdir}/pf_axi2dBW_HalfPlane_NR128.txt   silo redirect=tmp_
#
#
# First run some short simulations to make sure the early expansion is ok.
#
# N=064
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${ddir}/BWaxi2D_HalfPlane_n064_EINT_FKJav01 \
 redirect=${ddir}/msg_BWaxi2D_HalfPlane_n064_EINT_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
opfreq_time=1.0e20 finishtime=3.16e10 eqntype=9  &
#
# N=128
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${ddir}/BWaxi2D_HalfPlane_n128_EINT_FKJav01 \
 redirect=${ddir}/msg_BWaxi2D_HalfPlane_n128_EINT_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
opfreq_time=1.0e20 finishtime=3.16e10 eqntype=9  &
#
wait
#
# Then run simulations for the full 50kyr to make sure they all converge to roughly the same answer.
#
# N=064
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${ddir}/BWaxi2D_HalfPlane_n064_EINT_FKJav01 \
 redirect=${ddir}/msg_BWaxi2D_HalfPlane_n064_EINT_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
opfreq_time=1.0e20 eqntype=9   &
#
# N=128
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${ddir}/BWaxi2D_HalfPlane_n128_EINT_FKJav01 \
 redirect=${ddir}/msg_BWaxi2D_HalfPlane_n128_EINT_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
opfreq_time=1.0e20 eqntype=9   &
#
wait
# --------------------------------------------------------------------
# Now I should have a sequence for each simulation of:
# 0: initial conditions
# 1: Results at  1 kyr (3.16e10)
# 2: Results at 50 kyr (1.58e12)
# 3: checkpoint file 1 .9999998.silo
# 4: checkpoint file 2 .9999999.silo
#
# So let's get rid of the last two.
#
rm ${ddir}/BWaxi2D_HalfPlane*.999999*.silo
#
# Now we can make visit images, and files 1 and 2 will be the ones of
# interest.  I can tile them for each resolution.
# --------------------------------------------------------------------
# 
cd $tdir
#
infile_base=BWaxi2D_HalfPlane_n128_EINT_FKJav01
outfile_base=BWaxi2D_HalfPlane_n128_EINT_FKJav01
./make_BWaxi2D_image.sh $tdir $code_dir $ddir $visit_cmd $infile_base $outfile_base
#
infile_base=BWaxi2D_HalfPlane_n064_EINT_FKJav01
outfile_base=BWaxi2D_HalfPlane_n064_EINT_FKJav01
./make_BWaxi2D_image.sh $tdir $code_dir $ddir $visit_cmd $infile_base $outfile_base
#
convert ${ddir}/BWaxi2D_HalfPlane_n064_EINT_FKJav01_Dens_01.tif ${ddir}/BWaxi2D_HalfPlane_n064_EINT_FKJav01_Dens_02.tif +append EINT_n064_time.jpeg
convert ${ddir}/BWaxi2D_HalfPlane_n128_EINT_FKJav01_Dens_01.tif ${ddir}/BWaxi2D_HalfPlane_n128_EINT_FKJav01_Dens_02.tif +append EINT_n128_time.jpeg
convert EINT_n064_time.jpeg EINT_n128_time.jpeg -append EINT_res_time.jpeg
eog EINT_res_time.jpeg &


#./run_axisymmetric_BW.sh $tdir/blastwave_axi2d $code_dir ${ddir}/blastwave_axi2d
#./visit_BWaxi2d_make_plots.sh $tdir/blastwave_axi2d $code_dir ${ddir}/blastwave_axi2d $visit_cmd
#
# End of Test Problem: Blastwave (1D,2D (slab+axi),3D)
#---------------------------------------------------------------------
#---------------------------------------------------------------------

exit

#---------------------------------------------------------------------
# Run the Double Mach Reflection test
#---------------------------------------------------------------------
echo "EINT:: DMR TEST :: Starting"
cd $code_dir
make -f Makefile.serial.code; make -f Makefile.serial.icgenerator

#
# generate ICs
#
./icgen ${test_dir}/test_DMR/pf_DMR_n260.txt silo
./icgen ${test_dir}/test_DMR/pf_DMR_n520.txt silo
#
# Run low and high resolution models.
#
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_EINT_av10 cfl=0.5 \
    redirect=${data_dir}/msg_DMRm10t60_n260_EINT_av10 solver=3 eqntype=9 AVtype=1 EtaVisc=0.1 cfl=0.3 &
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_EINT_av10 cfl=0.5 \
    redirect=${data_dir}/msg_DMRm10t60_n520_EINT_av10 solver=3 eqntype=9 AVtype=1 EtaVisc=0.1 cfl=0.3 &
wait
#
# Generate figures to compare with other solvers
#
cd ${test_dir}/test_DMR/
./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_EINT_av10
./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_EINT_av10
cd ${data_dir}
## lo-res
convert -crop 955x325+175+225 DMRm10t60_n260_EINT_av10.tif DMRm10t60_n260_EINT_av10.tif 
convert DMRm10t60_n260_EINT_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, EINT/ETOT solver (RPV=3), A.V. eta=0.1' DMRm10t60_n260_EINT_av10.tif
convert DMRm10t60_n260_EINT_av10.tif DMRm10t60_n260_EINT_av10.jpeg
convert DMRm10t60_n260_EINT_av10.tif eps3:DMRm10t60_n260_EINT_av10.eps
## hi-res
convert -crop 955x325+175+225 DMRm10t60_n520_EINT_av10.tif DMRm10t60_n520_EINT_av10.tif 
convert DMRm10t60_n520_EINT_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, EINT/ETOT solver (RPV=3), A.V. eta=0.1' DMRm10t60_n520_EINT_av10.tif
convert DMRm10t60_n520_EINT_av10.tif DMRm10t60_n520_EINT_av10.jpeg
convert DMRm10t60_n520_EINT_av10.tif eps3:DMRm10t60_n520_EINT_av10.eps
#
# move figs to test-dir
#
mv DMR*.eps DMR*.jpeg ${test_dir}/test_DMR/
cd -
echo "DOUBLE MACH REFLECTION: FINISHED GENERATING FIGURES, OPENING VIEWER"
echo "COMPARE IMAGES TO: http://www.astro.uni-bonn.de/~jmackey/jmac/node10.html"
eog *.jpeg &
cd ${test_dir}
echo "EINT:: DMR TEST :: Done"
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# 1D/2D shock tubes
#---------------------------------------------------------------------
echo "************** RUNNING  1D/2D SHOCK-TUBES **************"
tdir=${test_dir}/test_ShockTubes
cd ${test_dir}/test_ShockTubes/
# WARNING: HI-RES SHOCK-TUBE TESTS TAKE HOURS TO FINISH.
#./run_hiresST_tests.sh    $test_dir $code_dir $data_dir/STtest
mkdir $data_dir/STtest
ddir=${data_dir}/STtest


cd ${code_dir}
# HYDRO
./icgen ${tdir}/pf_st_toro1_hires.txt silo redirect=${ddir}/ic_Toro1_
./icgen ${tdir}/pf_st_toro2_hires.txt silo redirect=${ddir}/ic_Toro2_
./icgen ${tdir}/pf_st_toro3_hires.txt silo redirect=${ddir}/ic_Toro3_
./icgen ${tdir}/pf_st_toro4_hires.txt silo redirect=${ddir}/ic_Toro4_
./icgen ${tdir}/pf_st_toro5_hires.txt silo redirect=${ddir}/ic_Toro5_
#
./icgen ${tdir}/pf_st_toro1.txt silo redirect=${ddir}/ic_Toro1_
./icgen ${tdir}/pf_st_toro2.txt silo redirect=${ddir}/ic_Toro2_
./icgen ${tdir}/pf_st_toro3.txt silo redirect=${ddir}/ic_Toro3_
./icgen ${tdir}/pf_st_toro4.txt silo redirect=${ddir}/ic_Toro4_
./icgen ${tdir}/pf_st_toro5.txt silo redirect=${ddir}/ic_Toro5_
#
# Now run the 1D shock-tube tests
#
./main_serial IC_Toro1.silo 5 1 outfile=${ddir}/Toro1_EINT opfreq=0 redirect=${ddir}/msg_Toro1_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3
./main_serial IC_Toro2.silo 5 1 outfile=${ddir}/Toro2_EINT opfreq=0 redirect=${ddir}/msg_Toro2_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3
./main_serial IC_Toro3.silo 5 1 outfile=${ddir}/Toro3_EINT opfreq=0 redirect=${ddir}/msg_Toro3_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3
./main_serial IC_Toro4.silo 5 1 outfile=${ddir}/Toro4_EINT opfreq=0 redirect=${ddir}/msg_Toro4_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3
./main_serial IC_Toro5.silo 5 1 outfile=${ddir}/Toro5_EINT opfreq=0 redirect=${ddir}/msg_Toro5_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3
#
./main_serial IC_Toro10k1.silo 5 1 outfile=${ddir}/Toro10k1_EINT opfreq=0 redirect=${ddir}/msg_Toro10k1_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3 &
./main_serial IC_Toro10k2.silo 5 1 outfile=${ddir}/Toro10k2_EINT opfreq=0 redirect=${ddir}/msg_Toro10k2_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3 &
./main_serial IC_Toro10k3.silo 5 1 outfile=${ddir}/Toro10k3_EINT opfreq=0 redirect=${ddir}/msg_Toro10k3_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3 &
./main_serial IC_Toro10k4.silo 5 1 outfile=${ddir}/Toro10k4_EINT opfreq=0 redirect=${ddir}/msg_Toro10k4_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3 &
./main_serial IC_Toro10k5.silo 5 1 outfile=${ddir}/Toro10k5_EINT opfreq=0 redirect=${ddir}/msg_Toro10k5_EINT_ artvisc=0.1 ooa=2 eqntype=9 solver=3 cfl=0.3 &
wait
#
# Now for 2D
#
tdir=${test_dir}/test_ShockTubes/2D
# HYDRO
./icgen ${tdir}/pf_st2D_toro1.txt silo redirect=${ddir}/ic2d_
./icgen ${tdir}/pf_st2D_toro2.txt silo redirect=${ddir}/ic2d_
./icgen ${tdir}/pf_st2D_toro3.txt silo redirect=${ddir}/ic2d_
./icgen ${tdir}/pf_st2D_toro4.txt silo redirect=${ddir}/ic2d_
./icgen ${tdir}/pf_st2D_toro5.txt silo redirect=${ddir}/ic2d_
#
# Now run the shock-tube tests
#
./main_serial IC_Toro2D1.silo 5 1 outfile=${ddir}/Toro2D1_EINT redirect=${ddir}/msg_Toro2D1_EINT artvisc=0.1 ooa=2 opfreq_time=5.0 eqntype=9 solver=3 cfl=0.3 &
./main_serial IC_Toro2D2.silo 5 1 outfile=${ddir}/Toro2D2_EINT redirect=${ddir}/msg_Toro2D2_EINT artvisc=0.1 ooa=2 opfreq_time=5.0 eqntype=9 solver=3 cfl=0.3 &
./main_serial IC_Toro2D3.silo 5 1 outfile=${ddir}/Toro2D3_EINT redirect=${ddir}/msg_Toro2D3_EINT artvisc=0.1 ooa=2 opfreq_time=5.0 eqntype=9 solver=3 cfl=0.3 &
./main_serial IC_Toro2D4.silo 5 1 outfile=${ddir}/Toro2D4_EINT redirect=${ddir}/msg_Toro2D4_EINT artvisc=0.1 ooa=2 opfreq_time=5.0 eqntype=9 solver=3 cfl=0.3 &
./main_serial IC_Toro2D5.silo 5 1 outfile=${ddir}/Toro2D5_EINT redirect=${ddir}/msg_Toro2D5_EINT artvisc=0.1 ooa=2 opfreq_time=5.0 eqntype=9 solver=3 cfl=0.3 &
wait
#
cd $tdir
make -f Makefile.twoD2oneD
file1=`ls ${ddir}/Toro2D1_EINT.0*[1-9]*.silo`
echo "./twoD2oneD_shocktubes $file1 silo -1 ${tdir}/Toro2D1_plotdata.txt 0.0 1.0 0.0 1.0"
./twoD2oneD_shocktubes $file1 silo -1 ${tdir}/Toro2D1_EINT_plotdata.txt 0.0 1.0 0.0 1.0
file1=`ls ${ddir}/Toro2D2_EINT.0*[1-9]*.silo`
./twoD2oneD_shocktubes $file1 silo -1 ${tdir}/Toro2D2_EINT_plotdata.txt 0.0 1.0 0.0 1.0
file1=`ls ${ddir}/Toro2D3_EINT.0*[1-9]*.silo`
./twoD2oneD_shocktubes $file1 silo -1 ${tdir}/Toro2D3_EINT_plotdata.txt 0.0 1.0 0.0 1.0
file1=`ls ${ddir}/Toro2D4_EINT.0*[1-9]*.silo`
./twoD2oneD_shocktubes $file1 silo -1 ${tdir}/Toro2D4_EINT_plotdata.txt 0.0 1.0 0.0 1.0
file1=`ls ${ddir}/Toro2D5_EINT.0*[1-9]*.silo`
./twoD2oneD_shocktubes $file1 silo -1 ${tdir}/Toro2D5_EINT_plotdata.txt 0.0 1.0 0.0 1.0
#

tdir=${test_dir}/test_ShockTubes
cd $tdir
#
# 1D plots
#
file1=`ls ${ddir}/Toro1_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/../../STtest/Toro10k1.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro1_EINT_TOT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro1_EINT_TOT

file1=`ls ${ddir}/Toro2_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/../../STtest/Toro10k2.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro2_EINT_TOT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro2_EINT_TOT

file1=`ls ${ddir}/Toro3_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/../../STtest/Toro10k3.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro3_EINT_TOT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro3_EINT_TOT

file1=`ls ${ddir}/Toro4_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/../../STtest/Toro10k4.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro4_EINT_TOT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro4_EINT_TOT

file1=`ls ${ddir}/Toro5_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/../../STtest/Toro10k5.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro5_EINT_TOT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro5_EINT_TOT
#
file1=`ls ${ddir}/Toro1_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/Toro10k1_EINT.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro1_EINT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro1_EINT

file1=`ls ${ddir}/Toro2_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/Toro10k2_EINT.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro2_EINT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro2_EINT

file1=`ls ${ddir}/Toro3_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/Toro10k3_EINT.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro3_EINT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro3_EINT

file1=`ls ${ddir}/Toro4_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/Toro10k4_EINT.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro4_EINT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro4_EINT

file1=`ls ${ddir}/Toro5_EINT.0*[1-9]*.txt`
file2=`ls ${ddir}/Toro10k5_EINT.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${tdir}/Toro5_EINT "
./plot_ST.sh $file1 $file2  ${tdir}/Toro5_EINT

#
# 2D plots
# Offsets are calculated from the interface at x0 via xoff = x0*(cos(t)-1) + 0.5*sin(t)
# Look in shock_tube.cc and the twoD2oneD.cc conversion program to see why.
#
file2=`ls ${ddir}/../../STtest/Toro10k1.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${tdir}/2D/Toro2D1_EINT_plotdata.txt ${tdir}/2D/Toro2D_T1_EINT 0.1919351
file2=`ls ${ddir}/../../STtest/Toro10k2.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${tdir}/2D/Toro2D2_EINT_plotdata.txt ${tdir}/2D/Toro2D_T2_EINT 0.1708204
file2=`ls ${ddir}/../../STtest/Toro10k3.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${tdir}/2D/Toro2D3_EINT_plotdata.txt ${tdir}/2D/Toro2D_T3_EINT 0.1708204
file2=`ls ${ddir}/../../STtest/Toro10k4.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${tdir}/2D/Toro2D4_EINT_plotdata.txt ${tdir}/2D/Toro2D_T4_EINT 0.1813777
file2=`ls ${ddir}/../../STtest/Toro10k5.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${tdir}/2D/Toro2D5_EINT_plotdata.txt ${tdir}/2D/Toro2D_T5_EINT 0.1391486
#
echo "FINISHED MAKING PLOTS.  OPENING A VIEWER NOW."
echo "COMPARE IMAGES TO http://www.astro.uni-bonn.de/~jmackey/jmac/node9.html"
eog *EINT*.jpeg &
eog 2D/*EINT*.jpeg &

#./run_ST2D_tests.sh $test_dir $code_dir $data_dir/STtest
#./run_ST_tests.sh   $test_dir $code_dir $data_dir/STtest
#./make_ST_plots.sh  ${test_dir}/test_ShockTubes $code_dir $data_dir/STtest
echo "************** FINISHED 1D/2D SHOCK-TUBES **************"
#---------------------------------------------------------------------

