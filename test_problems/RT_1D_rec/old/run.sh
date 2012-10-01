#!/bin/bash


# 2011.10.22 JM: Testing multifrequency photoionisation
../../bin/icgen_serial pf_rtt_Sph1D_ERG_n128_nh2_MP3_MultiFreq.txt silo
../../bin/icgen_serial pf_rtt_Sph1D_ERG_n128_nh2_MP4_MultiFreq.txt silo

../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2_MP3_MultiFreq.silo 5 1 outfile=./temp_outputs/mp3_multi_UV_cflp10 op_criterion=0 opfreq=1 finishtime=3.16e10 cfl=0.10 optype=silo
../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2_MP3_MultiFreq.silo 5 1 outfile=./temp_outputs/mp3_multi_UV_cflp01 op_criterion=0 opfreq=1 finishtime=3.16e10 cfl=0.01 optype=silo
#
../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2_MP4_MultiFreq.silo 5 1 outfile=./temp_outputs/mp4_multi_UV_cflp100 op_criterion=0 opfreq=1 finishtime=3.16e10 cfl=0.100 optype=silo
../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2_MP4_MultiFreq.silo 5 1 outfile=./temp_outputs/mp4_multi_UV_cflp010 op_criterion=0 opfreq=1 finishtime=3.16e10 cfl=0.010 optype=silo
../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2_MP4_MultiFreq.silo 5 1 outfile=./temp_outputs/mp4_multi_UV_cflp001 op_criterion=0 opfreq=1 finishtime=3.16e10 cfl=0.001 optype=silo
#
../RT_1D_norec/plot_radius temp_outputs/rad_mp3_multi_UV_cflp10 temp_outputs/mp3_multi_UV_cflp10 0 1 5 silo
../RT_1D_norec/plot_radius temp_outputs/rad_mp3_multi_UV_cflp01 temp_outputs/mp3_multi_UV_cflp01 0 1 5 silo
../RT_1D_norec/plot_radius temp_outputs/rad_mp4_multi_UV_cflp100 temp_outputs/mp4_multi_UV_cflp100 0 1 5 silo
../RT_1D_norec/plot_radius temp_outputs/rad_mp4_multi_UV_cflp010 temp_outputs/mp4_multi_UV_cflp010 0 1 5 silo
../RT_1D_norec/plot_radius temp_outputs/rad_mp4_multi_UV_cflp001 temp_outputs/mp4_multi_UV_cflp001 0 1 5 silo
#
exit

../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2_MP3_MultiFreq.silo 5 1 outfile=./temp_outputs/mpv3_MF_ns050v4 redirect=./temp_outputs/msg_mpv3_MF_ns050v4 op_criterion=1 opfreq_time=3.16e9 optype=6 checkpt_freq=100000
exit

#valgrind --tool=callgrind --dump-instr=yes ../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2.silo 5 1 outfile=./temp_outputs/mpv2 redirect=./temp_outputs/msg_mpv2 opfreq=5 optype=6 checkpt_freq=100000 finishtime=1.0e9

#valgrind --tool=callgrind --dump-instr=yes ../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2_MP3.silo 5 1 outfile=./temp_outputs/mpv3 redirect=./temp_outputs/msg_mpv3 opfreq=5 optype=6 checkpt_freq=100000 finishtime=1.0e9

#valgrind --tool=callgrind --dump-instr=yes ../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2_MP2_MultiFreq.silo 5 1 outfile=./temp_outputs/mpv2_MF redirect=./temp_outputs/msg_mpv2_MF opfreq=5 optype=6 checkpt_freq=100000 finishtime=1.0e9

valgrind --tool=callgrind --dump-instr=yes ../../bin/main_serial IC_rtt_Sph1D_ERG_n128_nh2_MP3_MultiFreq.silo 5 1 outfile=./temp_outputs/mpv3_MF redirect=./temp_outputs/msg_mpv3_MF opfreq=5 optype=6 checkpt_freq=100000 finishtime=1.0e9


