set terminal postscript enhanced color eps
#set terminal postscript enhanced eps
set size 0.7071,0.7071
set xlabel "Time (t_{rec})"
set key right bottom
unset log y
set ylabel "Ratio R_{a}/R_{if}(t)"
set yrange [0.95:1.02]
set log x
set xrange [0.01:30]
set title ""
set output "/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/RT_1D_rec/COMP1D_radius_error_nh2_n128_C2RAY.eps"
set label "1D, n_H=100 cm^{-3}, Alg1" at 0.015,1.0075
set label "Cell {/Symbol d}{/Symbol t}=30 for N=128 runs" at 0.015,1.015
g(x)=(1.0-exp(-x))**(1.0/3.0)
RH=6.0e18
RI=4.5172e18
tH=3.861e10
tI=3.861e10
f(x)=1.0
h(x)=1.0+4.6875e16/(RI*g(x))
j(x)=1.0-4.6875e16/(RI*g(x))
plot f(x) w l lw 2 lt -1 title 'Correct Solution', \
     h(x) w l lw 2 lt 0 title "One cell error", j(x) w l lw 2 lt 0 notitle, \
     './data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr030.txt' u ($1/tI):(($5/RI)/g($1/tI)) w lp lw 2 lt 1  title 'N=128, {/Symbol d}t= 0.30 t_{rec}', \
     './data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr010.txt' u ($1/tI):(($5/RI)/g($1/tI)) w lp lw 2 lt 2  title 'N=128, {/Symbol d}t= 0.10 t_{rec}', \
     './data/rtt_Sph1D_ERG_n128_nh2_C2ray_tr003.txt' u ($1/tI):(($5/RI)/g($1/tI)) w lp lw 2 lt 3  title 'N=128, {/Symbol d}t= 0.03 t_{rec}'
#pause -1
quit

