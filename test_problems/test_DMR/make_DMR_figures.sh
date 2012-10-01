#!/bin/bash

# call with ./make_DMR_figures.sh $test_dir $code_dir $data_dir
test_dir=$1
code_dir=$2
data_dir=$3
visit_cmd=$4

echo "GENERATING PLOTS"
cd ${test_dir}/test_DMR/

#./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_EINT_av10
#cd ${data_dir}
## lo-res
#convert -crop 955x325+175+225 DMRm10t60_n260_EINT_av10.tif DMRm10t60_n260_EINT_av10.tif 
#convert DMRm10t60_n260_EINT_av10.tif -pointsize 20 -annotate +250+275 \
# 'Nx=260, EINT/ETOT solver (RPV=3), A.V. eta=0.1' DMRm10t60_n260_EINT_av10.tif
#exit

./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_Hyb_av00 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_Hyb_av10 $visit_cmd

./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_RCV_av00 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_RCV_av10 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_RCV_Hcor $visit_cmd

./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_RPV_av00 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_RPV_av10 $visit_cmd

./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_FVS_av00 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n260_FVS_av10 $visit_cmd


#
# Crop figures down to the bare minimum of white space.
#
cd ${data_dir}
# lo-res
convert -crop 955x325+175+225 DMRm10t60_n260_Hyb_av00.tif DMRm10t60_n260_Hyb_av00.tif 
convert -crop 955x325+175+225 DMRm10t60_n260_Hyb_av10.tif DMRm10t60_n260_Hyb_av10.tif 

convert -crop 955x325+175+225 DMRm10t60_n260_RCV_av00.tif DMRm10t60_n260_RCV_av00.tif 
convert -crop 955x325+175+225 DMRm10t60_n260_RCV_av10.tif DMRm10t60_n260_RCV_av10.tif 
convert -crop 955x325+175+225 DMRm10t60_n260_RCV_Hcor.tif DMRm10t60_n260_RCV_Hcor.tif 

convert -crop 955x325+175+225 DMRm10t60_n260_RPV_av00.tif DMRm10t60_n260_RPV_av00.tif 
convert -crop 955x325+175+225 DMRm10t60_n260_RPV_av10.tif DMRm10t60_n260_RPV_av10.tif 

convert -crop 955x325+175+225 DMRm10t60_n260_FVS_av00.tif DMRm10t60_n260_FVS_av00.tif 
convert -crop 955x325+175+225 DMRm10t60_n260_FVS_av10.tif DMRm10t60_n260_FVS_av10.tif 
#
# Now annotate each figure with its attributes
#
convert DMRm10t60_n260_Hyb_av00.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, Primitive Var. Riemann solver, A.V. eta=0.0' DMRm10t60_n260_Hyb_av00.tif
convert DMRm10t60_n260_Hyb_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, Primitive Var. Riemann solver, A.V. eta=0.1' DMRm10t60_n260_Hyb_av10.tif

convert DMRm10t60_n260_RCV_av00.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, Roe Conserved Var. flux solver, A.V. eta=0.0' DMRm10t60_n260_RCV_av00.tif
convert DMRm10t60_n260_RCV_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, Roe Conserved Var. flux solver, A.V. eta=0.1' DMRm10t60_n260_RCV_av10.tif
convert DMRm10t60_n260_RCV_Hcor.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, Roe Conserved Var. flux solver, H-correction' DMRm10t60_n260_RCV_Hcor.tif

convert DMRm10t60_n260_RPV_av00.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, Roe Primitive Var. flux solver, A.V. eta=0.0' DMRm10t60_n260_RPV_av00.tif
convert DMRm10t60_n260_RPV_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, Roe Primitive Var. flux solver, A.V. eta=0.1' DMRm10t60_n260_RPV_av10.tif

convert DMRm10t60_n260_FVS_av00.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, Flux Vector Splitting   solver, A.V. eta=0.0' DMRm10t60_n260_FVS_av00.tif
convert DMRm10t60_n260_FVS_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=260, Flux Vector Splitting   solver, A.V. eta=0.1' DMRm10t60_n260_FVS_av10.tif

#
# Now convert to jpegs and eps.
#
# LO-RES
convert DMRm10t60_n260_Hyb_av00.tif DMRm10t60_n260_Hyb_av00.jpeg
convert DMRm10t60_n260_Hyb_av10.tif DMRm10t60_n260_Hyb_av10.jpeg
convert DMRm10t60_n260_RCV_av00.tif DMRm10t60_n260_RCV_av00.jpeg
convert DMRm10t60_n260_RCV_av10.tif DMRm10t60_n260_RCV_av10.jpeg
convert DMRm10t60_n260_RCV_Hcor.tif DMRm10t60_n260_RCV_Hcor.jpeg
convert DMRm10t60_n260_RPV_av00.tif DMRm10t60_n260_RPV_av00.jpeg
convert DMRm10t60_n260_RPV_av10.tif DMRm10t60_n260_RPV_av10.jpeg
convert DMRm10t60_n260_FVS_av00.tif DMRm10t60_n260_FVS_av00.jpeg
convert DMRm10t60_n260_FVS_av10.tif DMRm10t60_n260_FVS_av10.jpeg

convert DMRm10t60_n260_RCV_av00.tif eps2:DMRm10t60_n260_RCV_av00.eps
convert DMRm10t60_n260_RCV_av10.tif eps2:DMRm10t60_n260_RCV_av10.eps
convert DMRm10t60_n260_Hyb_av00.tif eps2:DMRm10t60_n260_Hyb_av00.eps
convert DMRm10t60_n260_Hyb_av10.tif eps2:DMRm10t60_n260_Hyb_av10.eps
convert DMRm10t60_n260_RCV_Hcor.tif eps2:DMRm10t60_n260_RCV_Hcor.eps 
convert DMRm10t60_n260_RPV_av00.tif eps3:DMRm10t60_n260_RPV_av00.eps
convert DMRm10t60_n260_RPV_av10.tif eps3:DMRm10t60_n260_RPV_av10.eps
convert DMRm10t60_n260_FVS_av00.tif eps3:DMRm10t60_n260_FVS_av00.eps
convert DMRm10t60_n260_FVS_av10.tif eps3:DMRm10t60_n260_FVS_av10.eps

#
# move figures to test_problems/test_DMR/
#
mv DMR*n260*.eps DMR*n260*.jpeg ${test_dir}/test_DMR/
cd -
echo "DOUBLE MACH REFLECTION: FINISHED GENERATING FIGURES, OPENING VIEWER"
echo "COMPARE IMAGES TO: http://www.astro.uni-bonn.de/~jmackey/jmac/node10.html"
eog *n260*.jpeg &
cd ${test_dir}

########## TEMP #########
#exit
########## TEMP #########


cd ${test_dir}/test_DMR/
./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_Hyb_av00 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_Hyb_av10 $visit_cmd

./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_RCV_av00 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_RCV_av10 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_RCV_Hcor $visit_cmd

./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_RPV_av00 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_RPV_av10 $visit_cmd

./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_FVS_av00 $visit_cmd
./plot_DMR_results.sh ${data_dir} DMRm10t60_n520_FVS_av10 $visit_cmd

cd ${data_dir}
# hi-res
convert -crop 955x325+175+225 DMRm10t60_n520_Hyb_av00.tif DMRm10t60_n520_Hyb_av00.tif 
convert -crop 955x325+175+225 DMRm10t60_n520_Hyb_av10.tif DMRm10t60_n520_Hyb_av10.tif 

convert -crop 955x325+175+225 DMRm10t60_n520_RCV_av00.tif DMRm10t60_n520_RCV_av00.tif 
convert -crop 955x325+175+225 DMRm10t60_n520_RCV_av10.tif DMRm10t60_n520_RCV_av10.tif 
convert -crop 955x325+175+225 DMRm10t60_n520_RCV_Hcor.tif DMRm10t60_n520_RCV_Hcor.tif 

convert -crop 955x325+175+225 DMRm10t60_n520_RPV_av00.tif DMRm10t60_n520_RPV_av00.tif 
convert -crop 955x325+175+225 DMRm10t60_n520_RPV_av10.tif DMRm10t60_n520_RPV_av10.tif 

convert -crop 955x325+175+225 DMRm10t60_n520_FVS_av00.tif DMRm10t60_n520_FVS_av00.tif 
convert -crop 955x325+175+225 DMRm10t60_n520_FVS_av10.tif DMRm10t60_n520_FVS_av10.tif 


convert DMRm10t60_n520_Hyb_av00.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, Primitive Var. Riemann solver, A.V. eta=0.0' DMRm10t60_n520_Hyb_av00.tif
convert DMRm10t60_n520_Hyb_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, Primitive Var. Riemann solver, A.V. eta=0.1' DMRm10t60_n520_Hyb_av10.tif

convert DMRm10t60_n520_RCV_av00.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, Roe Conserved Var. flux solver, A.V. eta=0.0' DMRm10t60_n520_RCV_av00.tif
convert DMRm10t60_n520_RCV_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, Roe Conserved Var. flux solver, A.V. eta=0.1' DMRm10t60_n520_RCV_av10.tif
convert DMRm10t60_n520_RCV_Hcor.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, Roe Conserved Var. flux solver, H-correction' DMRm10t60_n520_RCV_Hcor.tif

convert DMRm10t60_n520_RPV_av00.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, Roe Primitive Var. flux solver, A.V. eta=0.0' DMRm10t60_n520_RPV_av00.tif
convert DMRm10t60_n520_RPV_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, Roe Primitive Var. flux solver, A.V. eta=0.1' DMRm10t60_n520_RPV_av10.tif

convert DMRm10t60_n520_FVS_av00.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, Flux Vector Splitting   solver, A.V. eta=0.0' DMRm10t60_n520_FVS_av00.tif
convert DMRm10t60_n520_FVS_av10.tif -pointsize 20 -annotate +250+275 \
 'Nx=520, Flux Vector Splitting   solver, A.V. eta=0.1' DMRm10t60_n520_FVS_av10.tif

# HI-RES
convert DMRm10t60_n520_Hyb_av00.tif DMRm10t60_n520_Hyb_av00.jpeg
convert DMRm10t60_n520_Hyb_av10.tif DMRm10t60_n520_Hyb_av10.jpeg

convert DMRm10t60_n520_RCV_av00.tif DMRm10t60_n520_RCV_av00.jpeg
convert DMRm10t60_n520_RCV_av10.tif DMRm10t60_n520_RCV_av10.jpeg
convert DMRm10t60_n520_RCV_Hcor.tif DMRm10t60_n520_RCV_Hcor.jpeg

convert DMRm10t60_n520_RPV_av00.tif DMRm10t60_n520_RPV_av00.jpeg
convert DMRm10t60_n520_RPV_av10.tif DMRm10t60_n520_RPV_av10.jpeg

convert DMRm10t60_n520_FVS_av00.tif DMRm10t60_n520_FVS_av00.jpeg
convert DMRm10t60_n520_FVS_av10.tif DMRm10t60_n520_FVS_av10.jpeg

#
convert DMRm10t60_n520_RCV_av00.tif eps2:DMRm10t60_n520_RCV_av00.eps
convert DMRm10t60_n520_RCV_av10.tif eps2:DMRm10t60_n520_RCV_av10.eps
convert DMRm10t60_n520_Hyb_av00.tif eps2:DMRm10t60_n520_Hyb_av00.eps
convert DMRm10t60_n520_Hyb_av10.tif eps2:DMRm10t60_n520_Hyb_av10.eps
convert DMRm10t60_n520_RCV_Hcor.tif eps2:DMRm10t60_n520_RCV_Hcor.eps 

convert DMRm10t60_n520_RPV_av00.tif eps3:DMRm10t60_n520_RPV_av00.eps
convert DMRm10t60_n520_RPV_av10.tif eps3:DMRm10t60_n520_RPV_av10.eps
convert DMRm10t60_n520_FVS_av00.tif eps3:DMRm10t60_n520_FVS_av00.eps
convert DMRm10t60_n520_FVS_av10.tif eps3:DMRm10t60_n520_FVS_av10.eps


#
# move figures to test_problems/test_DMR/
#
mv DMR*n520*.eps DMR*n520*.jpeg ${test_dir}/test_DMR/
cd -
echo "DOUBLE MACH REFLECTION: FINISHED GENERATING FIGURES, OPENING VIEWER"
echo "COMPARE IMAGES TO: http://www.astro.uni-bonn.de/~jmackey/jmac/node10.html"
eog *n520*.jpeg &
cd ${test_dir}/test_DMR

mv DMR*.jpeg DMR*.eps ${test_dir}/

exit

