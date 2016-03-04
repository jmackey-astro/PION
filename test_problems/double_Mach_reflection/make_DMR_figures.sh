#!/bin/bash

# call with ./make_DMR_figures.sh $test_dir $code_dir $data_dir $visit_cmd $resolution
test_dir=$1
code_dir=$2
data_dir=$3
visit_cmd=$4
RESOLUTION=$5

##########  FOR A GIVEN RESOLUTION  ##########
echo "GENERATING PLOTS for resolution ${RESOLUTION}"
cd ${test_dir}

./plot_DMR_results.sh ${data_dir} DMRm10t60_n${RESOLUTION}_Hyb_av10 $visit_cmd
sleep 10
./plot_DMR_results.sh ${data_dir} DMRm10t60_n${RESOLUTION}_RCV_av10 $visit_cmd
sleep 10
./plot_DMR_results.sh ${data_dir} DMRm10t60_n${RESOLUTION}_RCV_Hcor $visit_cmd
sleep 10
./plot_DMR_results.sh ${data_dir} DMRm10t60_n${RESOLUTION}_RPV_av10 $visit_cmd
sleep 10
./plot_DMR_results.sh ${data_dir} DMRm10t60_n${RESOLUTION}_FVS_av10 $visit_cmd
sleep 10


#
# Crop figures down to the bare minimum of white space.
#
cd ${data_dir}
# lo-res
convert -crop 955x325+175+225 DMRm10t60_n${RESOLUTION}_Hyb_av10.png DMRm10t60_n${RESOLUTION}_Hyb_av10.png 
convert -crop 955x325+175+225 DMRm10t60_n${RESOLUTION}_RCV_av10.png DMRm10t60_n${RESOLUTION}_RCV_av10.png 
convert -crop 955x325+175+225 DMRm10t60_n${RESOLUTION}_RCV_Hcor.png DMRm10t60_n${RESOLUTION}_RCV_Hcor.png 
convert -crop 955x325+175+225 DMRm10t60_n${RESOLUTION}_RPV_av10.png DMRm10t60_n${RESOLUTION}_RPV_av10.png 
convert -crop 955x325+175+225 DMRm10t60_n${RESOLUTION}_FVS_av10.png DMRm10t60_n${RESOLUTION}_FVS_av10.png 
#
# Now annotate each figure with its attributes
#
convert DMRm10t60_n${RESOLUTION}_Hyb_av10.png -pointsize 20 -annotate +250+275 \
 "Nx=${RESOLUTION}, Primitive Var. Riemann solver, A.V. eta=0.1" DMRm10t60_n${RESOLUTION}_Hyb_av10.png
convert DMRm10t60_n${RESOLUTION}_RCV_av10.png -pointsize 20 -annotate +250+275 \
 "Nx=${RESOLUTION}, Roe Conserved Var. flux solver, A.V. eta=0.1" DMRm10t60_n${RESOLUTION}_RCV_av10.png
convert DMRm10t60_n${RESOLUTION}_RCV_Hcor.png -pointsize 20 -annotate +250+275 \
 "Nx=${RESOLUTION}, Roe Conserved Var. flux solver, H-correction" DMRm10t60_n${RESOLUTION}_RCV_Hcor.png
convert DMRm10t60_n${RESOLUTION}_RPV_av10.png -pointsize 20 -annotate +250+275 \
 "Nx=${RESOLUTION}, Roe Primitive Var. flux solver, A.V. eta=0.1" DMRm10t60_n${RESOLUTION}_RPV_av10.png
convert DMRm10t60_n${RESOLUTION}_FVS_av10.png -pointsize 20 -annotate +250+275 \
 "Nx=${RESOLUTION}, Flux Vector Splitting   solver, A.V. eta=0.1" DMRm10t60_n${RESOLUTION}_FVS_av10.png

mkdir -p ../FIGS
mv *.png ../FIGS/

exit
##########  FOR A GIVEN RESOLUTION  ##########

