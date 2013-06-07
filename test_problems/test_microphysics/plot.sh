#!/bin/bash

cat << EOF > temp.gp
set terminal postscript enhanced color eps font "Times,20" \
 fontfile "/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb"
#set size 0.7071,0.7071
set size 1.0,1.0
set border linewidth 1.5
set tics scale 1.5
set lmargin 8.0
set tmargin 0.5
set rmargin 1.25
set bmargin 3.5
set key at 1400,7e-21
#set key top left
set log y
set log x
set xrange [10:25000]
set yrange [1e-28:1e-20]
set xlabel "T (K)" offset 0,0
set ylabel "Heating/cooling [erg cm^{-3} s^{-1}]" offset 2,0
set output "fig_cool_MPv3.eps"
plot   "./MPv3_heating_cooling_nH1_r10.txt" u 1:(abs(\$7)) w l lt 1 lc 0 lw 4 title "H II, n_H=1 cm^{-3},   f(H+)=0.9964", \
       "./MPv3_heating_cooling_nH1_r10.txt" u 1:7          w p pt 1 lc 0 ps 1.5 notitle, \
       "./MPv3_heating_cooling_nH0.1_r30.txt" u 1:(abs(\$2)) w l lt 5 lc 0 lw 3 title "WNM, n_H=0.1 cm^{-3}, f(H+)=2.1e-2", \
       "./MPv3_heating_cooling_nH0.1_r30.txt" u 1:2          w p pt 1 lc 0 ps 1.5 notitle, \
       "./MPv3_heating_cooling_nH1_r10.txt" u 1:(abs(\$2)) w l lt 5 lc 3 lw 3 title "WNM, n_H=1 cm^{-3},   f(H+)=4.0e-3", \
       "./MPv3_heating_cooling_nH1_r10.txt" u 1:2          w p pt 1 lc 3 ps 1.5 notitle, \
       "./MPv3_heating_cooling_nH10_r3.txt" u 1:(abs(\$2)) w l lt 5 lc 2 lw 3 title "WNM, n_H=10 cm^{-3},  f(H+)=5.2e-4", \
       "./MPv3_heating_cooling_nH10_r3.txt" u 1:(\$2)          w p pt 1 lc 2 ps 1.5 notitle, \
       "./MPv3_heating_cooling_nH100_r1.txt" u 1:(abs(\$2)) w l lt 5 lc 5 lw 3 title "WNM, n_H=100 cm^{-3}, f(H+)=7.5e-5", \
       "./MPv3_heating_cooling_nH100_r1.txt" u 1:(\$2)          w p pt 1 lc 5 ps 1.5 notitle
exit
EOF
gnuplot temp.gp


exit

plot   "./MPv3_heating_cooling_nH1_r10.txt" u 1:(abs(\$7)) w l lt 1 lc 0 lw 4 title "Unshielded, f(H+)=0.9964", \
       "./MPv3_heating_cooling_nH1_r10.txt" u 1:7          w p pt 1 lc 0 ps 1.5 notitle, \
       "./MPv3_heating_cooling_nH1_r10.txt" u 1:(abs(\$2)) w l lt 5 lc 0 lw 3 title "Shielded, n_H=1,   f(H+)=4.0e-3", \
       "./MPv3_heating_cooling_nH1_r10.txt" u 1:2          w p pt 1 lc 0 ps 1.5 notitle, \
       "./MPv3_heating_cooling_nH10_r3.txt" u 1:(abs(\$2)*1.0e-2) w l lt 5 lc 2 lw 3 title "Shielded, n_H=10,  f(H+)=5.2e-4", \
       "./MPv3_heating_cooling_nH10_r3.txt" u 1:(\$2*1.0e-2)          w p pt 1 lc 2 ps 1.5 notitle, \
       "./MPv3_heating_cooling_nH100_r1.txt" u 1:(abs(\$2)*1.0e-4) w l lt 5 lc 3 lw 3 title "Shielded, n_H=100, f(H+)=7.5e-5", \
       "./MPv3_heating_cooling_nH100_r1.txt" u 1:(\$2*1.0e-4)          w p pt 1 lc 3 ps 1.5 notitle
exit

## old plot!
cat << EOF > temp.gp
set terminal postscript enhanced color eps font "Times,20" \
 fontfile "/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb"
set size 1.0,1.0
#set size 0.7071,0.7071
set lmargin 8.0
set tmargin 0.5
set rmargin 1.25
set bmargin 3.5
set key at 400,2e-28
set log y
set log x
set xrange [10:25000]
set yrange [1e-30:1e-19]
set xlabel "T (K)" offset 0,0
set ylabel "Heating/cooling [erg.cm^{-3}.s^{-1}]" offset 2,0
set output "fig_cool_MPv3.eps"
plot   "./MPv3heating_cooling.txt" u 1:(abs(\$7)) w l lt 1 lc 0 lw 2 title "Unshielded, f(H+)=0.995", \
       "./MPv3heating_cooling.txt" u 1:7          w p pt 1 lc 0 notitle, \
       "./MPv3heating_cooling.txt" u 1:(abs(\$6)) w l lt 1 lc 2 lw 2 title "Unshielded, f(H+)=0.100", \
       "./MPv3heating_cooling.txt" u 1:6          w p pt 1 lc 2 notitle, \
       "./MPv3heating_cooling.txt" u 1:(abs(\$4)) w l lt 5 lc 0 lw 3 title "Shielded,   f(H+)=0.995", \
       "./MPv3heating_cooling.txt" u 1:4          w p pt 1 lc 0 notitle, \
       "./MPv3heating_cooling.txt" u 1:(abs(\$3)) w l lt 5 lc 2 lw 3 title "Shielded,   f(H+)=0.100", \
       "./MPv3heating_cooling.txt" u 1:3          w p pt 1 lc 2 notitle, \
       "./MPv3heating_cooling.txt" u 1:(abs(\$2)) w l lt 5 lc 3 lw 3 title "Shielded,   f(H+)=0.002", \
       "./MPv3heating_cooling.txt" u 1:2          w p pt 1 lc 3 notitle
exit
EOF
gnuplot temp.gp


exit


       "./MPv3heating_cooling.txt" u 1:(abs(\$5)) w l lt 1 lc 3 lw 2 title "Unshielded, f(H+)=0.002", \
       "./MPv3heating_cooling.txt" u 1:5          w p pt 1 lc 3 notitle, \

