#!/bin/bash

# generate ics
../RT_1D_norec/icgen_serial pf_rtt_Sph1D_ERG_n128_nh2_C2ray.txt silo

# run models with different timestep limits.
../RT_1D_norec/main_serial_C2RAY_030tr IC_rtt_Sph1D_ERG_n128_nh2_C2ray.silo 5 1 \
 outfile=./data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr030 optype=6 op_criterion=0 opfreq=1
../RT_1D_norec/main_serial_C2RAY_010tr IC_rtt_Sph1D_ERG_n128_nh2_C2ray.silo 5 1 \
 outfile=./data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr010 optype=6 op_criterion=0 opfreq=1
../RT_1D_norec/main_serial_C2RAY_003tr IC_rtt_Sph1D_ERG_n128_nh2_C2ray.silo 5 1 \
 outfile=./data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr003 optype=6 op_criterion=0 opfreq=1

# calculate IF radius as function of time.
../RT_1D_norec/plot_radius ./data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr030 \
   ./data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr030 0 1 5 silo
../RT_1D_norec/plot_radius ./data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr010 \
   ./data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr010 0 1 5 silo
../RT_1D_norec/plot_radius ./data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr003 \
   ./data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr003 0 1 5 silo

# make plot
gnuplot make_c2ray_fig.gnuplot
convert -density 300 -quality 100 COMP1D_radius_error_nh2_n128_C2RAY.eps \
                                  COMP1D_radius_error_nh2_n128_C2RAY.jpeg

# delete all the data files
rm ./data/rtt_Sph1D_ERG_n128_nh2_C2ray*

exit
