#!/bin/bash
#

g++ -Wall -g -DSERIAL parallelionisation.cc ../../testing/global.cc ../../testing/uniformGrid.cc ../../testing/dataio.cc -lreadline -lcfitsio
./a.out d1_const_parallel/nh0ss1e7_n1k_dt1_acc3sv2 ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt1_acc3sv2 0 1 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt2_acc3sv2 ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt2_acc3sv2 0 1 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt3_acc3sv2 ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt3_acc3sv2 0 8 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt4_acc3sv2 ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt4_acc3sv2 0 50 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt1_acc3sv2 ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt1_acc3sv2 0 1 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt2_acc3sv2 ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt2_acc3sv2 0 1 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt3_acc3sv2 ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt3_acc3sv2 0 8 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt4_acc3sv2 ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt4_acc3sv2 0 50 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt1_acc3sv2 ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt1_acc3sv2 0 1 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt2_acc3sv2 ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt2_acc3sv2 0 1 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt3_acc3sv2 ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt3_acc3sv2 0 8 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt4_acc3sv2 ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt4_acc3sv2 0 50 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt1_acc3sv2 ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt1_acc3sv2 0 1 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt2_acc3sv2 ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt2_acc3sv2 0 1 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt3_acc3sv2 ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt3_acc3sv2 0 8 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt4_acc3sv2 ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt4_acc3sv2 0 50 5
exit


############################################################
# Monday 26/5/08                                           #
# Analysis for 1D photo-ionisation problems, w/recombs.    #
############################################################
g++ -Wall -g -DSERIAL parallelionisation.cc ../../testing/global.cc ../../testing/uniformGrid.cc ../../testing/dataio.cc -lreadline -lcfitsio
#############
## nh=1
#  t/dt=10
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt1_acc3hard   ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt1_acc3h 0 1 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt1_acc3simple ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt1_acc3s 0 1 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt1_acc4hard   ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt1_acc4h 0 1 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt1_acc4simple ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt1_acc4s 0 1 5
#  t/dt=10^2
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt2_acc3hard   ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt2_acc3h 0 1 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt2_acc3simple ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt2_acc3s 0 1 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt2_acc4hard   ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt2_acc4h 0 1 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt2_acc4simple ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt2_acc4s 0 1 5
#  t/dt=10^3
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt3_acc3hard   ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt3_acc3h 0 8 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt3_acc3simple ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt3_acc3s 0 8 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt3_acc4hard   ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt3_acc4h 0 8 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt3_acc4simple ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt3_acc4s 0 8 5
#  t/dt=10^4
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt4_acc3hard   ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt4_acc3h 0 50 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt4_acc3simple ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt4_acc3s 0 50 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt4_acc4hard   ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt4_acc4h 0 50 5
./a.out d1_const_parallel/RRnh0ss2e6_n1c_dt4_acc4simple ../../results/PE1dcart_parallel_RRnh0ss2e6_n1c_dt4_acc4s 0 50 5
#############
## nh=10
#  t/dt=10
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt1_acc3hard   ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt1_acc3h 0 1 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt1_acc3simple ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt1_acc3s 0 1 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt1_acc4hard   ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt1_acc4h 0 1 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt1_acc4simple ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt1_acc4s 0 1 5
#  t/dt=10^2
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt2_acc3hard   ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt2_acc3h 0 1 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt2_acc3simple ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt2_acc3s 0 1 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt2_acc4hard   ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt2_acc4h 0 1 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt2_acc4simple ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt2_acc4s 0 1 5
#  t/dt=10^3
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt3_acc3hard   ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt3_acc3h 0 8 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt3_acc3simple ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt3_acc3s 0 8 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt3_acc4hard   ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt3_acc4h 0 8 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt3_acc4simple ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt3_acc4s 0 8 5
#  t/dt=10^4
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt4_acc3hard   ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt4_acc3h 0 50 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt4_acc3simple ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt4_acc3s 0 50 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt4_acc4hard   ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt4_acc4h 0 50 5
./a.out d1_const_parallel/RRnh1ss2e8_n1c_dt4_acc4simple ../../results/PE1dcart_parallel_RRnh1ss2e8_n1c_dt4_acc4s 0 50 5
#############
## nh=100
#  t/dt=10
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt1_acc3hard   ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt1_acc3h 0 1 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt1_acc3simple ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt1_acc3s 0 1 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt1_acc4hard   ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt1_acc4h 0 1 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt1_acc4simple ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt1_acc4s 0 1 5
#  t/dt=10^2
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt2_acc3hard   ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt2_acc3h 0 1 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt2_acc3simple ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt2_acc3s 0 1 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt2_acc4hard   ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt2_acc4h 0 1 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt2_acc4simple ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt2_acc4s 0 1 5
#  t/dt=10^3
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt3_acc3hard   ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt3_acc3h 0 8 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt3_acc3simple ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt3_acc3s 0 8 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt3_acc4hard   ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt3_acc4h 0 8 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt3_acc4simple ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt3_acc4s 0 8 5
#  t/dt=10^4
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt4_acc3hard   ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt4_acc3h 0 50 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt4_acc3simple ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt4_acc3s 0 50 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt4_acc4hard   ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt4_acc4h 0 50 5
./a.out d1_const_parallel/RRnh2ss2e10_n1c_dt4_acc4simple ../../results/PE1dcart_parallel_RRnh2ss2e10_n1c_dt4_acc4s 0 50 5
#############
## nh=1000
#  t/dt=10
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt1_acc3hard   ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt1_acc3h 0 1 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt1_acc3simple ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt1_acc3s 0 1 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt1_acc4hard   ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt1_acc4h 0 1 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt1_acc4simple ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt1_acc4s 0 1 5
#  t/dt=10^2
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt2_acc3hard   ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt2_acc3h 0 1 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt2_acc3simple ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt2_acc3s 0 1 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt2_acc4hard   ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt2_acc4h 0 1 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt2_acc4simple ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt2_acc4s 0 1 5
#  t/dt=10^3
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt3_acc3hard   ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt3_acc3h 0 8 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt3_acc3simple ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt3_acc3s 0 8 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt3_acc4hard   ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt3_acc4h 0 8 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt3_acc4simple ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt3_acc4s 0 8 5
#  t/dt=10^4
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt4_acc3hard   ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt4_acc3h 0 50 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt4_acc3simple ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt4_acc3s 0 50 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt4_acc4hard   ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt4_acc4h 0 50 5
./a.out d1_const_parallel/RRnh3ss2e12_n1c_dt4_acc4simple ../../results/PE1dcart_parallel_RRnh3ss2e12_n1c_dt4_acc4s 0 50 5
exit


############################################################
# Monday 19/5/08                                           #
# Analysis for 1D photo-ionisation problems, no recombs.   #
############################################################
g++ -Wall -g -DSERIAL parallelionisation.cc ../../testing/global.cc ../../testing/uniformGrid.cc ../../testing/dataio.cc -lreadline -lcfitsio
#############
## nh=1
#  t/dt=10
./a.out d1_const_parallel/nh0ss1e7_n1k_dt1_acc3hard   ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt1_acc3h 0 1 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt1_acc3simple ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt1_acc3s 0 1 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt1_acc4hard   ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt1_acc4h 0 1 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt1_acc4simple ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt1_acc4s 0 1 5
#  t/dt=10^2
./a.out d1_const_parallel/nh0ss1e7_n1k_dt2_acc3hard   ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt2_acc3h 0 1 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt2_acc3simple ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt2_acc3s 0 1 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt2_acc4hard   ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt2_acc4h 0 1 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt2_acc4simple ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt2_acc4s 0 1 5
#  t/dt=10^3
./a.out d1_const_parallel/nh0ss1e7_n1k_dt3_acc3hard   ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt3_acc3h 0 8 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt3_acc3simple ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt3_acc3s 0 8 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt3_acc4hard   ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt3_acc4h 0 8 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt3_acc4simple ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt3_acc4s 0 8 5
#  t/dt=10^4
./a.out d1_const_parallel/nh0ss1e7_n1k_dt4_acc3hard   ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt4_acc3h 0 50 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt4_acc3simple ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt4_acc3s 0 50 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt4_acc4hard   ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt4_acc4h 0 50 5
./a.out d1_const_parallel/nh0ss1e7_n1k_dt4_acc4simple ../../results/PE1dcart_parallel_nh0ss1e7_n1k_dt4_acc4s 0 50 5
#############
## nh=10
#  t/dt=10
./a.out d1_const_parallel/nh1ss1e7_n1k_dt1_acc3hard   ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt1_acc3h 0 1 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt1_acc3simple ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt1_acc3s 0 1 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt1_acc4hard   ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt1_acc4h 0 1 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt1_acc4simple ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt1_acc4s 0 1 5
#  t/dt=10^2
./a.out d1_const_parallel/nh1ss1e7_n1k_dt2_acc3hard   ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt2_acc3h 0 1 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt2_acc3simple ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt2_acc3s 0 1 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt2_acc4hard   ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt2_acc4h 0 1 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt2_acc4simple ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt2_acc4s 0 1 5
#  t/dt=10^3
./a.out d1_const_parallel/nh1ss1e7_n1k_dt3_acc3hard   ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt3_acc3h 0 8 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt3_acc3simple ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt3_acc3s 0 8 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt3_acc4hard   ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt3_acc4h 0 8 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt3_acc4simple ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt3_acc4s 0 8 5
#  t/dt=10^4
./a.out d1_const_parallel/nh1ss1e7_n1k_dt4_acc3hard   ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt4_acc3h 0 50 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt4_acc3simple ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt4_acc3s 0 50 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt4_acc4hard   ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt4_acc4h 0 50 5
./a.out d1_const_parallel/nh1ss1e7_n1k_dt4_acc4simple ../../results/PE1dcart_parallel_nh1ss1e7_n1k_dt4_acc4s 0 50 5
#############
## nh=100
#  t/dt=10
./a.out d1_const_parallel/nh2ss1e7_n1k_dt1_acc3hard   ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt1_acc3h 0 1 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt1_acc3simple ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt1_acc3s 0 1 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt1_acc4hard   ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt1_acc4h 0 1 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt1_acc4simple ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt1_acc4s 0 1 5
#  t/dt=10^2
./a.out d1_const_parallel/nh2ss1e7_n1k_dt2_acc3hard   ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt2_acc3h 0 1 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt2_acc3simple ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt2_acc3s 0 1 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt2_acc4hard   ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt2_acc4h 0 1 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt2_acc4simple ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt2_acc4s 0 1 5
#  t/dt=10^3
./a.out d1_const_parallel/nh2ss1e7_n1k_dt3_acc3hard   ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt3_acc3h 0 8 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt3_acc3simple ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt3_acc3s 0 8 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt3_acc4hard   ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt3_acc4h 0 8 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt3_acc4simple ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt3_acc4s 0 8 5
#  t/dt=10^4
./a.out d1_const_parallel/nh2ss1e7_n1k_dt4_acc3hard   ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt4_acc3h 0 50 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt4_acc3simple ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt4_acc3s 0 50 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt4_acc4hard   ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt4_acc4h 0 50 5
./a.out d1_const_parallel/nh2ss1e7_n1k_dt4_acc4simple ../../results/PE1dcart_parallel_nh2ss1e7_n1k_dt4_acc4s 0 50 5
#############
## nh=1000
#  t/dt=10
./a.out d1_const_parallel/nh3ss1e7_n1k_dt1_acc3hard   ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt1_acc3h 0 1 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt1_acc3simple ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt1_acc3s 0 1 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt1_acc4hard   ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt1_acc4h 0 1 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt1_acc4simple ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt1_acc4s 0 1 5
#  t/dt=10^2
./a.out d1_const_parallel/nh3ss1e7_n1k_dt2_acc3hard   ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt2_acc3h 0 1 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt2_acc3simple ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt2_acc3s 0 1 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt2_acc4hard   ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt2_acc4h 0 1 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt2_acc4simple ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt2_acc4s 0 1 5
#  t/dt=10^3
./a.out d1_const_parallel/nh3ss1e7_n1k_dt3_acc3hard   ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt3_acc3h 0 8 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt3_acc3simple ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt3_acc3s 0 8 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt3_acc4hard   ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt3_acc4h 0 8 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt3_acc4simple ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt3_acc4s 0 8 5
#  t/dt=10^4
./a.out d1_const_parallel/nh3ss1e7_n1k_dt4_acc3hard   ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt4_acc3h 0 50 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt4_acc3simple ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt4_acc3s 0 50 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt4_acc4hard   ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt4_acc4h 0 50 5
./a.out d1_const_parallel/nh3ss1e7_n1k_dt4_acc4simple ../../results/PE1dcart_parallel_nh3ss1e7_n1k_dt4_acc4s 0 50 5
exit

