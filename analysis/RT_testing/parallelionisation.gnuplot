
######################################################
## 27/5/08 Plots with just acc3s, with all timesteps.
######################################################
## all densities, multiple t/dt per plot, all at 'acc3simple'
set terminal postscript enhanced color eps
set xlabel "Time (Seconds)"
set ylabel "Radius Error (dx, abs.value)"
set key bottom right
set grid
set log y
set yrange [1e-2:10]
set xrange [0:2e13]
#set title "Nx=100, nh=1 per cc"
set output "d1_const_parallel/raderr_RRnh0ss2e06_nx100_dt.eps"
plot "d1_const_parallel/RRnh0ss2e6_n1c_dt1_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt2_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt3_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt4_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^4"
#pause -1
set yrange [1e-2:10]
set xrange [0:2e12]
#set title "Nx=100, nh=10 per cc"
set output "d1_const_parallel/raderr_RRnh1ss2e08_nx100_dt.eps"
plot "d1_const_parallel/RRnh1ss2e8_n1c_dt1_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt2_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt3_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt4_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^4"
#pause -1
set xrange [0:2e11]
#set title "Nx=100, nh=100 per cc"
set output "d1_const_parallel/raderr_RRnh2ss2e10_nx100_dt.eps"
plot "d1_const_parallel/RRnh2ss2e10_n1c_dt1_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt2_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt3_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt4_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^4"
#pause -1
set xrange [0:2e10]
#set title "Nx=100, nh=1000 per cc"
set output "d1_const_parallel/raderr_RRnh3ss2e12_nx100_dt.eps"
plot "d1_const_parallel/RRnh3ss2e12_n1c_dt1_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt2_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt3_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt4_acc3simple.txt" u 1:(abs($4-$3)/1.e17) w l  title "acc3s, t/dt=10^4"
#pause -1
set yrange [1e-2:1]
set xrange [0:2e10]
#set title "Nx=20, nh=1000 per cc"
set output "d1_const_parallel/raderr_RRnh3ss2e12_nx020_dt.eps"
plot "d1_const_parallel/RRnh3ss2e12_n20_dt1_acc3simple.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/RRnh3ss2e12_n20_dt2_acc3simple.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/RRnh3ss2e12_n20_dt3_acc3simple.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/RRnh3ss2e12_n20_dt4_acc3simple.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s, t/dt=10^4"
#pause -1
set yrange [1e-2:1]
set xrange [0:2e12]
#set title "Nx=20, nh=10 per cc"
set output "d1_const_parallel/raderr_RRnh1ss2e08_nx020_dt.eps"
plot "d1_const_parallel/RRnh1ss2e8_n20_dt1_acc3simple.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/RRnh1ss2e8_n20_dt2_acc3simple.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/RRnh1ss2e8_n20_dt3_acc3simple.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/RRnh1ss2e8_n20_dt4_acc3simple.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s, t/dt=10^4"
#pause -1
set xrange [0:2e13]
#set title "Nx=20, nh=10 per cc"
set output "d1_const_parallel/raderr_RRnh1ss2e08_nx020_dtlong.eps"
plot "d1_const_parallel/RRnh1ss2e8_n20_dt1_acc3simplev2.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s2, t/dt=10^1", \
     "d1_const_parallel/RRnh1ss2e8_n20_dt2_acc3simplev2.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s2, t/dt=10^2", \
     "d1_const_parallel/RRnh1ss2e8_n20_dt3_acc3simplev2.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s2, t/dt=10^3", \
     "d1_const_parallel/RRnh1ss2e8_n20_dt4_acc3simplev2.txt" u 1:(abs($4-$3)/5.e17) w l  title "acc3s2, t/dt=10^4"
#pause -1

###############################
######## NO RECOMBS!!! SAME PLOTS AS ABOVE
## all densities, multiple t/dt per plot, all at 'acc3simple'
set yrange [*:*]
set xrange [*:*]
set key top left
set xrange [0:1.e12]
set yrange [1e-3:20]
#set title "NO RECOMBS, Nx=1000, nh=1.587 per cc, Tau=0.1 per cell"
set output "d1_const_parallel/raderr_nh0ss1e7_nx1k_dt_old.eps"
plot "d1_const_parallel/nh0ss1e7_n1k_dt1_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/nh0ss1e7_n1k_dt2_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/nh0ss1e7_n1k_dt3_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/nh0ss1e7_n1k_dt4_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^4"
#pause -1
set yrange [*:*]
set yrange [1e-3:20]
#set title "NO RECOMBS, Nx=1000, nh=1.587 per cc, Tau=0.1 per cell, V2!!!"
set output "d1_const_parallel/raderr_nh0ss1e7_nx1k_dt.eps"
plot "d1_const_parallel/nh0ss1e7_n1k_dt1_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^1", \
     "d1_const_parallel/nh0ss1e7_n1k_dt2_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^2", \
     "d1_const_parallel/nh0ss1e7_n1k_dt3_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^3", \
     "d1_const_parallel/nh0ss1e7_n1k_dt4_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^4"
#pause -1
set yrange [*:*]
set xrange [0:1e13]
set yrange [-1:20]
set yrange [1e-4:20]
#set title "NO RECOMBS, Nx=1000, nh=15.87 per cc, Tau=1 per cell"
set output "d1_const_parallel/raderr_nh1ss1e7_nx1k_dt_old.eps"
plot "d1_const_parallel/nh1ss1e7_n1k_dt1_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/nh1ss1e7_n1k_dt2_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/nh1ss1e7_n1k_dt3_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/nh1ss1e7_n1k_dt4_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^4"
#pause -1
set yrange [*:*]
set yrange [1e-4:20]
#set title "NO RECOMBS, Nx=1000, nh=15.87 per cc, Tau=1 per cell  V2!!!!!!"
set output "d1_const_parallel/raderr_nh1ss1e7_nx1k_dt.eps"
plot "d1_const_parallel/nh1ss1e7_n1k_dt1_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^1", \
     "d1_const_parallel/nh1ss1e7_n1k_dt2_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^2", \
     "d1_const_parallel/nh1ss1e7_n1k_dt3_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^3", \
     "d1_const_parallel/nh1ss1e7_n1k_dt4_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^4"
#pause -1
set xrange [0:1e14]
set yrange [-1:20]
#set yrange [*:*]
set yrange [1e-4:50]
#set title "NO RECOMBS, Nx=1000, nh=158.7 per cc, Tau=10 per cell"
set output "d1_const_parallel/raderr_nh2ss1e7_nx1k_dt_old.eps"
plot "d1_const_parallel/nh2ss1e7_n1k_dt1_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/nh2ss1e7_n1k_dt2_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/nh2ss1e7_n1k_dt3_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/nh2ss1e7_n1k_dt4_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^4"
#pause -1
set yrange [-20:5]
set yrange [1e-4:50]
#set title "NO RECOMBS, Nx=1000, nh=158.7 per cc, Tau=10 per cell, V2!!!!"
set output "d1_const_parallel/raderr_nh2ss1e7_nx1k_dt.eps"
plot "d1_const_parallel/nh2ss1e7_n1k_dt1_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^1", \
     "d1_const_parallel/nh2ss1e7_n1k_dt2_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^2", \
     "d1_const_parallel/nh2ss1e7_n1k_dt3_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^3", \
     "d1_const_parallel/nh2ss1e7_n1k_dt4_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^4"
#pause -1
# the above plot is with a different approximation for the implicit int(exp(-tau),dt) step, using 0.5(t0+t1)hh
# instead of t1*hh, so it should converge faster, and it does for t/dt=100, but t/dt is still bad, only the
# other way now.  but this is good -- i'll use it.  it's not much more work, and i get an extra order of 
# accuracy in terms of convergence, which I need!
set xrange [0:1e15]
set yrange [-1:20]
set yrange [1e-4:50]
#set title "NO RECOMBS, Nx=1000, nh=1587 per cc, Tau=100 per cell"
set output "d1_const_parallel/raderr_nh3ss1e7_nx1k_dt_old.eps"
plot "d1_const_parallel/nh3ss1e7_n1k_dt1_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^1", \
     "d1_const_parallel/nh3ss1e7_n1k_dt2_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^2", \
     "d1_const_parallel/nh3ss1e7_n1k_dt3_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^3", \
     "d1_const_parallel/nh3ss1e7_n1k_dt4_acc3simple.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s, t/dt=10^4"
#pause -1
set yrange [*:*]
set yrange [1e-4:50]
#set title "NO RECOMBS, Nx=1000, nh=1587 per cc, Tau=100 per cell, V2!!!"
set output "d1_const_parallel/raderr_nh3ss1e7_nx1k_dt.eps"
plot "d1_const_parallel/nh3ss1e7_n1k_dt1_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^1", \
     "d1_const_parallel/nh3ss1e7_n1k_dt2_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^2", \
     "d1_const_parallel/nh3ss1e7_n1k_dt3_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^3", \
     "d1_const_parallel/nh3ss1e7_n1k_dt4_acc3sv2.txt" u 1:(abs($4-$3)/1.e16) w l  title "acc3s2, t/dt=10^4"
#pause -1
quit


#############################
## 22/5/08 1D plots w/RRec ##
#############################
set xlabel "Time (Seconds)"
set ylabel "Radius Error (cells)"
set key bottom right
set grid
############################
## nh=1
set yrange [-3:7]
#set yrange [*:*]
set title "Nx=100, nh=1 per cc, t/dt=10"
plot "d1_const_parallel/RRnh0ss2e6_n1c_dt1_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt1_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt1_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt1_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=1 per cc, t/dt=100"
plot "d1_const_parallel/RRnh0ss2e6_n1c_dt2_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt2_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt2_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt2_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=1 per cc, t/dt=1000"
plot "d1_const_parallel/RRnh0ss2e6_n1c_dt3_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt3_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt3_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt3_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=1 per cc, t/dt=10000"
plot "d1_const_parallel/RRnh0ss2e6_n1c_dt4_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt4_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt4_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh0ss2e6_n1c_dt4_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
############################
## nh=10
#set log y
set yrange [-5:1]
#set yrange [*:*]
set title "Nx=100, nh=10 per cc, t/dt=10"
plot "d1_const_parallel/RRnh1ss2e8_n1c_dt1_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt1_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt1_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt1_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=10 per cc, t/dt=100"
plot "d1_const_parallel/RRnh1ss2e8_n1c_dt2_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt2_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt2_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt2_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=10 per cc, t/dt=1000"
plot "d1_const_parallel/RRnh1ss2e8_n1c_dt3_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt3_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt3_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt3_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=10 per cc, t/dt=10000"
plot "d1_const_parallel/RRnh1ss2e8_n1c_dt4_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt4_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt4_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh1ss2e8_n1c_dt4_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
############################
## nh=100
set yrange [-6:1]
#set yrange [*:*]
set title "Nx=100, nh=100 per cc, t/dt=10"
plot "d1_const_parallel/RRnh2ss2e10_n1c_dt1_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt1_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt1_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt1_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=100 per cc, t/dt=100"
plot "d1_const_parallel/RRnh2ss2e10_n1c_dt2_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt2_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt2_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt2_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=100 per cc, t/dt=1000"
plot "d1_const_parallel/RRnh2ss2e10_n1c_dt3_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt3_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt3_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt3_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=100 per cc, t/dt=10000"
plot "d1_const_parallel/RRnh2ss2e10_n1c_dt4_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt4_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt4_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh2ss2e10_n1c_dt4_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
############################
## nh=1000
set yrange [-6:1]
#set yrange [*:*]
set title "Nx=100, nh=1000 per cc, t/dt=10"
plot "d1_const_parallel/RRnh3ss2e12_n1c_dt1_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt1_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt1_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt1_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=1000 per cc, t/dt=100"
plot "d1_const_parallel/RRnh3ss2e12_n1c_dt2_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt2_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt2_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt2_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=1000 per cc, t/dt=1000"
plot "d1_const_parallel/RRnh3ss2e12_n1c_dt3_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt3_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt3_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt3_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1
set title "Nx=100, nh=1000 per cc, t/dt=10000"
plot "d1_const_parallel/RRnh3ss2e12_n1c_dt4_acc3simple.txt" u 1:($4-$3)/1.e17 w l  title "acc3, simple int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt4_acc3hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc3, better int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt4_acc4simple.txt" u 1:($4-$3)/1.e17 w l  title "acc4, simple int.", \
     "d1_const_parallel/RRnh3ss2e12_n1c_dt4_acc4hard.txt"   u 1:($4-$3)/1.e17 w lp title "acc4, better int."
pause -1

#quit




################################
## 22/5/08 1D plots  NO RECOMBINATIONS!!! ##
##############################
set xlabel "Time (Seconds)"
set ylabel "Radius Error (cells)"
set key top left
set grid
############################
## t/dt=10
set yrange [-1:30]
set title "n1000, nh=1.587 per cc, t/dt=10"
plot "d1_const_parallel/nh0ss1e7_n1k_dt1_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt1_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt1_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt1_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=15.87 per cc, t/dt=10"
plot "d1_const_parallel/nh1ss1e7_n1k_dt1_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt1_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt1_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt1_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=158.7 per cc, t/dt=10"
plot "d1_const_parallel/nh2ss1e7_n1k_dt1_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt1_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt1_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt1_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=1587 per cc, t/dt=10"
plot "d1_const_parallel/nh3ss1e7_n1k_dt1_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt1_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt1_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt1_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
############################
## t/dt=100
#set log y
set yrange [-0.2:15]
set title "n1000, nh=1.587 per cc, t/dt=100"
plot "d1_const_parallel/nh0ss1e7_n1k_dt2_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt2_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt2_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt2_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=15.87 per cc, t/dt=100"
plot "d1_const_parallel/nh1ss1e7_n1k_dt2_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt2_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt2_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt2_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=158.7 per cc, t/dt=100"
plot "d1_const_parallel/nh2ss1e7_n1k_dt2_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt2_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt2_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt2_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=1587 per cc, t/dt=100"
plot "d1_const_parallel/nh3ss1e7_n1k_dt2_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt2_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt2_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt2_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
############################
## t/dt=10^3
set yrange [-0.2:5]
set title "n1000, nh=1.587 per cc, t/dt=1000"
plot "d1_const_parallel/nh0ss1e7_n1k_dt3_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt3_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt3_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt3_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=15.87 per cc, t/dt=1000"
plot "d1_const_parallel/nh1ss1e7_n1k_dt3_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt3_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt3_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt3_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=158.7 per cc, t/dt=1000"
plot "d1_const_parallel/nh2ss1e7_n1k_dt3_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt3_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt3_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt3_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=1587 per cc, t/dt=1000"
plot "d1_const_parallel/nh3ss1e7_n1k_dt3_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt3_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt3_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt3_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
############################
## t/dt=10^4
set yrange [-0.2:3.5]
set title "n1000, nh=1.587 per cc, t/dt=10000"
plot "d1_const_parallel/nh0ss1e7_n1k_dt4_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt4_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt4_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh0ss1e7_n1k_dt4_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=15.87 per cc, t/dt=10000"
plot "d1_const_parallel/nh1ss1e7_n1k_dt4_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt4_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt4_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh1ss1e7_n1k_dt4_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=158.7 per cc, t/dt=10000"
plot "d1_const_parallel/nh2ss1e7_n1k_dt4_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt4_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt4_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh2ss1e7_n1k_dt4_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
set title "n1000, nh=1587 per cc, t/dt=10000"
plot "d1_const_parallel/nh3ss1e7_n1k_dt4_acc3simple.txt" u 1:($4-$3)/1.e16 w l  title "acc3, simple int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt4_acc3hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc3, better int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt4_acc4simple.txt" u 1:($4-$3)/1.e16 w l  title "acc4, simple int.", \
     "d1_const_parallel/nh3ss1e7_n1k_dt4_acc4hard.txt"   u 1:($4-$3)/1.e16 w lp title "acc4, better int."
pause -1
quit
# hard integration is using the substepping for exp(-tau) in the analytic
# integration.  simple is just taking the initial value and hoping for the 
# best.  
# For these problems, simple is more accurate -- has no bias at t/dt=10^4,
# and less bias (40%) at t/dt=10^2
# so should just check it for a denser sim.
# dt4 has no bias with simple integration at err=3 and 4, whereas dt2 is out by
# 2 cells at err=4 and 3.5 cells at err=3.
# for dt1, however, the hard integration has less bias, by a lot for acc4 (not
# much improvement at acc3.
# What does all this mean???  At very low time resolution (IF crosses 100 cells per 
# step) the hard integration is better at high accuracy, but makes little difference
# at acc3.  At higher time resolution, the hard integration is actually less 
# accurate (but always monotonically decreasing error with increasing time resolution).
# So is a 14cells out of 1000 error acceptable for when we are crossing 100cells/step???
# If so, then easy integration is better at higher time resolution, and I should use it.

################################
## 21/5/08 2D plots ##
##############################
set xlabel "Time (Seconds)"
set ylabel "Radius Error (cells)"

set title "n500, nh=10.7 per cc"
plot "d2_const_parallel/ss1e11_n500_dt020.txt" u 1:(($100-$101)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt100.txt" u 1:(($100-$101)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt200.txt" u 1:(($100-$101)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt500.txt" u 1:(($100-$101)/2.e17) w lp
pause -1
set title "n500, nh=9.7 per cc"
plot "d2_const_parallel/ss1e11_n500_dt020.txt" u 1:(($98-$99)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt100.txt" u 1:(($98-$99)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt200.txt" u 1:(($98-$99)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt500.txt" u 1:(($98-$99)/2.e17) w lp
pause -1
set title "n500, nh=1.08 per cc"
plot "d2_const_parallel/ss1e11_n500_dt020.txt" u 1:(($50-$51)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt100.txt" u 1:(($50-$51)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt200.txt" u 1:(($50-$51)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt500.txt" u 1:(($50-$51)/2.e17) w lp
pause -1
set title "n250, nh=10.7 per cc"
plot "d2_const_parallel/ss1e11_n250_dt020.txt" u 1:(($100-$101)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n250_dt100.txt" u 1:(($100-$101)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n250_dt200.txt" u 1:(($100-$101)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n250_dt500.txt" u 1:(($100-$101)/4.e17) w lp
pause -1
set title "v2 files use better implicit method for exp(-tau), and it shows!"
plot "d2_const_parallel/ss1e11_n250_dt200.txt" u 1:(($100-$101)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n250_dt200_v2.txt" u 1:(($100-$101)/4.e17) w lp
pause -1
set title "n250, nh=1.08 per cc"
plot "d2_const_parallel/ss1e11_n250_dt020.txt" u 1:(($50-$51)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n250_dt100.txt" u 1:(($50-$51)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n250_dt200.txt" u 1:(($50-$51)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n250_dt500.txt" u 1:(($50-$51)/4.e17) w lp
pause -1
set title "n=various, nh=10.7 per cc"
plot "d2_const_parallel/ss1e11_n250_dt100.txt"  u 1:(($100-$101)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt100.txt"  u 1:(($100-$101)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n1000_dt100.txt" u 1:(($100-$101)/1.e17) w lp
pause -1
set title "n=various, nh=10.7 per cc"
set xlabel "Radius Error (cm)"
plot "d2_const_parallel/ss1e11_n250_dt100.txt"  u 1:(($100-$101)) w lp, \
     "d2_const_parallel/ss1e11_n500_dt100.txt"  u 1:(($100-$101)) w lp, \
     "d2_const_parallel/ss1e11_n1000_dt100.txt" u 1:(($100-$101)) w lp
pause -1
plot "d2_const_parallel/ss1e11_n250_dt500.txt" u 1:(($100-$101)) w lp, \
     "d2_const_parallel/ss1e11_n500_dt500.txt" u 1:(($100-$101)) w lp, \
     "d2_const_parallel/ss1e11_n1000_dt500.txt" u 1:(($100-$101)) w lp
pause -1
set title "v2 files use better implicit method for exp(-tau), and it shows!"
set key top left
plot "d2_const_parallel/ss1e11_n250_dt200.txt"  u 1:(($100-$101)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n250_dt200_v2.txt"  u 1:(($100-$101)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt200.txt"  u 1:(($100-$101)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt200_v2.txt"  u 1:(($100-$101)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n1000_dt200_v2.txt" u 1:(($100-$101)/1.e17) w lp
pause -1
set title "v2 files use better implicit method for exp(-tau), and it shows!"
set key top left
plot "d2_const_parallel/ss1e11_n250_dt200.txt"  u 1:(($50-$51)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n250_dt200_v2.txt"  u 1:(($50-$51)/4.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt200.txt"  u 1:(($50-$51)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n500_dt200_v2.txt"  u 1:(($50-$51)/2.e17) w lp, \
     "d2_const_parallel/ss1e11_n1000_dt200_v2.txt" u 1:(($50-$51)/1.e17) w lp
pause -1
