#!/bin/bash

DDIR=/vol/jraid/jraid_1/jmackey/StarBench/DType1D
#DDIR=/vol/klaipeda3/scratch/jmackey/Spitzer_siloV2

rm fig/*.*

TTn01024=`ls ${DDIR}/txt_DTypeEarly_D1_TTI_n001024.0*.txt | head -n21 | tail -n1`
HCn01024=`ls ${DDIR}/txt_DTypeEarly_D1_HC_n001024.0*.txt | head -n21 | tail -n1`

echo "**** T=0.12 Myr****"
grep time $TTn01024
grep time $HCn01024

cat << EOF > temp.gp
#set terminal postscript enhanced color eps dashlength 2.0 font "Times-Roman,22" \
# fontfile "/usr/share/texlive/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy10.pfb"
set terminal postscript enhanced color eps font "Arial,26" \
  fontfile "/usr/share/texlive/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy10.pfb"

# sample line for controlling the key font and line sample.
#set key samplen 2 spacing .5 font ",8"
set key spacing 1.3

set size 1.0,1.0  # 0.7071,0.7071
set size 1,1
set ylabel "Shell Mass (M_{/CMSY10 \014})" offset 3,0
set log xy
set grid

set style line 11 lt 1 lc rgb "red" lw 4 pt 1 ps 1.5
set style line 12 lt 2 lc rgb "red" lw 4 pt 2 ps 1.5
set style line 13 lt 3 lc rgb "red" lw 4 pt 3 ps 1.5
set style line 14 lt 4 lc rgb "red" lw 4 pt 4 ps 1.5
set style line 15 lt 5 lc rgb "red" lw 4 pt 5 ps 1.5
set style line 16 lt 6 lc rgb "red" lw 4 pt 6 ps 1.5

set style line 21 lt 1 lc rgb "gray20" lw 3 pt 1 ps 1.5
set style line 22 lt 2 lc rgb "gray20" lw 3 pt 2 ps 1.5
set style line 23 lt 3 lc rgb "gray20" lw 3 pt 3 ps 1.5
set style line 24 lt 4 lc rgb "gray20" lw 3 pt 4 ps 1.5
set style line 25 lt 5 lc rgb "gray20" lw 3 pt 5 ps 1.5
set style line 26 lt 6 lc rgb "gray20" lw 3 pt 6 ps 1.5

set style line 31 lt 1 lc rgb "blue" lw 4 pt 1 ps 1.5
set style line 32 lt 2 lc rgb "blue" lw 4 pt 2 ps 1.5
set style line 33 lt 3 lc rgb "blue" lw 4 pt 3 ps 1.5
set style line 34 lt 4 lc rgb "blue" lw 4 pt 4 ps 1.5
set style line 35 lt 5 lc rgb "blue" lw 4 pt 5 ps 1.5
set style line 36 lt 6 lc rgb "blue" lw 4 pt 6 ps 1.5

set style line 41 lt 1 lc rgb "web-green" lw 4 pt 1 ps 1.5
set style line 42 lt 2 lc rgb "web-green" lw 4 pt 2 ps 1.5
set style line 43 lt 3 lc rgb "web-green" lw 4 pt 3 ps 1.5
set style line 44 lt 4 lc rgb "web-green" lw 4 pt 4 ps 1.5
set style line 45 lt 5 lc rgb "web-green" lw 4 pt 5 ps 1.5
set style line 46 lt 6 lc rgb "web-green" lw 4 pt 6 ps 1.5

set border lw 2

set size 1.0,1.0
#set size 0.7071,0.7071
set tmargin 1.5
set rmargin 1.25
set bmargin 2.5
set key top left 
set yrange [*:*]
unset log y
unset log x
set xlabel "r (pc)" offset 0,0.6


set output "fig/DTE1D_HCTT_Density_t0p1Myr.eps"
set xrange [1.08:1.155]
set yrange [5e1:2e5]
set ylabel "n_{H} (cm^{-3})" offset 3.5,0
set ylabel "" offset 3.5,0
#set title "Density, Temperature at t = 0.1 Myr" offset 0,-0.6
set log y
set lmargin 7.5
set key top left 
unset grid
plot "${TTn01024}" u (\$1):(\$2/1.67e-24) w lp ls 31 title "n_{H} (cm^{-3}) TTI", \
     "${TTn01024}" u (\$1):(\$4) w l ls 11 title "T(K) TTI", \
     "${HCn01024}" u (\$1):(\$2/1.67e-24) w lp ls 33 title "n_{H} (cm^{-3}) NEQ", \
     "${HCn01024}" u (\$1):(\$4) w l ls 13 title "T(K) NEQ"


EOF


echo "**** T=0.05 Myr****"
TTn01024=`ls ${DDIR}/txt_DTypeEarly_D1_TTI_n001024.0*.txt | head -n11 | tail -n1`
HCn01024=`ls ${DDIR}/txt_DTypeEarly_D1_HC_n001024.0*.txt | head -n11 | tail -n1`

echo "**** T=0.05 Myr****"
grep time $TTn01024
grep time $HCn01024


cat << EOF >> temp.gp

set output "fig/DTE1D_HCTT_Density_t0p05Myr.eps"
set xrange [0.74:0.85]
set yrange [5e1:2e5]
set ylabel "n_{H} (cm^{-3})" offset 3.5,0
set ylabel "" offset 3.5,0
#set title "Density, Temperature at t = 0.05 Myr" offset 0,-0.6
set log y
set lmargin 7.5
set key top left 
unset grid
plot "${TTn01024}" u (\$1):(\$2/1.67e-24) w lp ls 31 title "n_{H} (cm^{-3}) TTI", \
     "${TTn01024}" u (\$1):(\$4) w l ls 11 title "T(K) TTI", \
     "${HCn01024}" u (\$1):(\$2/1.67e-24) w lp ls 33 title "n_{H} (cm^{-3}) NEQ", \
     "${HCn01024}" u (\$1):(\$4) w l ls 13 title "T(K) NEQ"
EOF

echo "**** T=0.025 Myr****"
TTn01024=`ls ${DDIR}/txt_DTypeEarly_D1_TTI_n001024.0*.txt | head -n6 | tail -n1`
HCn01024=`ls ${DDIR}/txt_DTypeEarly_D1_HC_n001024.0*.txt | head -n6 | tail -n1`

echo "**** T=0.025 Myr****"
grep time $TTn01024
grep time $HCn01024


cat << EOF >> temp.gp

set output "fig/DTE1D_HCTT_Density_t0p025Myr.eps"
set xrange [0.74:0.85]
set xrange [0.535:0.615]
set yrange [5e1:2e5]
set ylabel "n_{H} (cm^{-3})" offset 3.5,0
set ylabel "" offset 3.5,0
#set title "Density, Temperature at t = 0.025 Myr" offset 0,-0.6
set log y
set lmargin 7.5
set key top left 
unset grid
plot "${TTn01024}" u (\$1):(\$2/1.67e-24) w lp ls 31 title "n_{H} (cm^{-3}) TTI", \
     "${TTn01024}" u (\$1):(\$4) w l ls 11 title "T(K) TTI", \
     "${HCn01024}" u (\$1):(\$2/1.67e-24) w lp ls 33 title "n_{H} (cm^{-3}) NEQ", \
     "${HCn01024}" u (\$1):(\$4) w l ls 13 title "T(K) NEQ"

#plot "${TTn01024}" u (\$1):(\$2/1.67e-24) w lp ls 31 title "n_{H} (cm^{-3}) TTI", \
#     "${TTn01024}" u (\$1):(\$4) w l ls 11 title "T(K) TTI", \
#     "${TTn01024}" u (\$1):(\$5*1.0e5) w l ls 21 title "x(H^{+})x10^5 TTI", \
#     "${HCn01024}" u (\$1):(\$2/1.67e-24) w lp ls 32 title "n_{H} (cm^{-3}) NEQ", \
#     "${HCn01024}" u (\$1):(\$4) w l ls 12 title "T(K) NEQ", \
#     "${HCn01024}" u (\$1):(\$5*1.0e5) w l ls 22 title "x(H^{+})x10^5 NEQ"


quit
EOF

gnuplot temp.gp
~/active/bin/eps2jpeg.sh fig/DTE
exit


