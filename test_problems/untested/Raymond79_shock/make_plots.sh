#!/bin/bash


#f1=RSH1D_n128_v100_Ray79E_HHe_C4.0023052.txt
#f2=RSH1D_n256_v100_Ray79E_HHe_C4.0046090.txt
#f3=RSH1D_n512_v100_Ray79E_HHe_C4.0091419.txt
f1=RSH1D_n128_v100_Ray79E_HHe_C4.00022800.txt
f2=RSH1D_n256_v100_Ray79E_HHe_C4.00045799.txt
f3=RSH1D_n512_v100_Ray79E_HHe_C4.00091293.txt
f4=RSH1D_n1024_v100_Ray79E_HHe_C4.00182641.txt
offset=31.35

cat << EOF  > gnu2.plt
set terminal postscript enhanced color eps
set size 1.0,0.5
set xrange [24:34]
set yrange [-0.1:1.5]
set xlabel "position (10^{15} cm)" 0.0,0.5
set ylabel "ion fraction" 1.0,0.0

set output "n128_He.eps"
plot '${f1}' u (\$1/1.0e15):12 w l title "He(1+)", \
 '${f1}' u (\$1/1.0e15):13 w l title "He(2+)", \
 '${f1}' u (\$1/1.0e15):(1-\$12-\$13) w l title "He(0)", \
 '${f1}' u (\$1/1.0e15):(\$14/1.0e5) w l lt -1 title "T(10^5K)"
set output "n128_H.eps"
plot '${f1}' u (\$1/1.0e15):10 w l title "x(e-)", \
     '${f1}' u (\$1/1.0e15):11 w l title "x(H+)", \
'${f1}' u (\$1/1.0e15):(1-\$11) w l title "x(H0)", \
'${f1}' u (\$1/1.0e15):(\$14/1.0e5) w l lt -1 title "T(10^5K)"

set output "n256_He.eps"
plot '${f2}' u (\$1/1.0e15):12 w l title "He(1+)", \
 '${f2}' u (\$1/1.0e15):13 w l title "He(2+)", \
 '${f2}' u (\$1/1.0e15):(1-\$12-\$13) w l title "He(0)", \
 '${f2}' u (\$1/1.0e15):(\$14/1.0e5) w l lt -1 title "T(10^5K)"
set output "n256_H.eps"
plot '${f2}' u (\$1/1.0e15):10 w l title "x(e-)", \
     '${f2}' u (\$1/1.0e15):11 w l title "x(H+)", \
'${f2}' u (\$1/1.0e15):(1-\$11) w l title "x(H0)", \
'${f2}' u (\$1/1.0e15):(\$14/1.0e5) w l lt -1 title "T(10^5K)"

set output "n512_He.eps"
plot '${f3}' u (\$1/1.0e15):12 w l title "He(1+)", \
 '${f3}' u (\$1/1.0e15):13 w l title "He(2+)", \
 '${f3}' u (\$1/1.0e15):(1-\$12-\$13) w l title "He(0)", \
 '${f3}' u (\$1/1.0e15):(\$14/1.0e5) w l lt -1 title "T(10^5K)"
set output "n512_H.eps"
plot '${f3}' u (\$1/1.0e15):10 w l title "x(e-)", \
     '${f3}' u (\$1/1.0e15):11 w l title "x(H+)", \
     '${f3}' u (\$1/1.0e15):(1-\$11) w l title "x(H0)", \
     '${f3}' u (\$1/1.0e15):(\$14/1.0e5) w l lt -1 title "T(10^5K)"

set output "n1024_He.eps"
#set xrange [31:32]
plot '${f4}' u (\$1/1.0e15):12 w l title "He(1+)", \
 '${f4}' u (\$1/1.0e15):13 w l title "He(2+)", \
 '${f4}' u (\$1/1.0e15):(1-\$12-\$13) w l title "He(0)", \
 '${f4}' u (\$1/1.0e15):(\$14/1.0e5) w l lt -1 title "T(10^5K)"
set output "n1024_H.eps"
plot '${f4}' u (\$1/1.0e15):10 w l title "x(e-)", \
     '${f4}' u (\$1/1.0e15):11 w l title "x(H+)", \
     '${f4}' u (\$1/1.0e15):(1-\$11) w l title "x(H0)", \
     '${f4}' u (\$1/1.0e15):(\$14/1.0e5) w l lt -1 title "T(10^5K)"

set xrange [-2:10]
set yrange [-0.1:1.75]
set output "H_res12.eps"
plot '${f2}' u (${offset}-\$1/1.0e15-0.25):10 w l lt 1 lw 2 title "x(e-)", \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):11 w l lt 2 lw 2 title "x(H+)", \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):(1-\$11) w l lt 3 lw 2 title "x(H0)", \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):(\$14/1.0e5) w l lt 0 lw 2 title "T(10^5K)", \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):10 w p lt 1 pt 1 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):11 w p lt 2 pt 2 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):(1-\$11) w p lt 3 pt 3 notitle
#, \
#     '${f1}' u (${offset}-\$1/1.0e15-0.75):(\$14/1.0e5) w lp lt 4 pt 4 notitle

set output "H_res23.eps"
plot '${f3}' u (${offset}-\$1/1.0e15):10 w l lt 1 lw 2 title "x(e-)", \
     '${f3}' u (${offset}-\$1/1.0e15):11 w l lt 2 lw 2 title "x(H+)", \
     '${f3}' u (${offset}-\$1/1.0e15):(1-\$11) w l lt 3 lw 2 title "x(H0)", \
     '${f3}' u (${offset}-\$1/1.0e15):(\$14/1.0e5) w l lt 0 lw 2 title "T(10^5K)", \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):10 w p lt 1 pt 1 notitle, \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):11 w p lt 2 pt 2 notitle, \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):(1-\$11) w p lt 3 pt 3 notitle
#, \
#     '${f2}' u (${offset}-\$1/1.0e15-0.25):(\$14/1.0e5) w lp lt 4 pt 4 notitle

set output "H_res14.eps"
plot '${f4}' u (${offset}-\$1/1.0e15+0.1):10 w l lt 1 lw 2 title "x(e-)", \
     '${f4}' u (${offset}-\$1/1.0e15+0.1):11 w l lt 2 lw 2 title "x(H+)", \
     '${f4}' u (${offset}-\$1/1.0e15+0.1):(1-\$11) w l lt 3 lw 2 title "x(H0)", \
     '${f4}' u (${offset}-\$1/1.0e15+0.1):(\$14/1.0e5) w l lt 0 lw 2 title "T(10^5K)", \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):10 w p lt 1 pt 1 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):11 w p lt 2 pt 2 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):(1-\$11) w p lt 3 pt 3 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):(\$14/1.0e5) w l lt 0 lw 2 lc 4 notitle

set output "H_res13.eps"
plot '${f3}' u (${offset}-\$1/1.0e15):10 w l lt 1 lw 2 title "x(e-)", \
     '${f3}' u (${offset}-\$1/1.0e15):11 w l lt 2 lw 2 title "x(H+)", \
     '${f3}' u (${offset}-\$1/1.0e15):(1-\$11) w l lt 3 lw 2 title "x(H0)", \
     '${f3}' u (${offset}-\$1/1.0e15):(\$14/1.0e5) w l lt 0 lw 2 title "T(10^5K)", \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):10 w p lt 1 pt 1 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):11 w p lt 2 pt 2 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):(1-\$11) w p lt 3 pt 3 notitle
#, \
#     '${f1}' u (${offset}-\$1/1.0e15-0.75):(\$14/1.0e5) w lp lt 4 pt 4 notitle

set output "He_res12.eps"
plot '${f2}' u (${offset}-\$1/1.0e15-0.25):12 w l lt 1 lw 2 title "He(1+)", \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):13 w l lt 2 lw 2 title "He(2+)", \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):(1-\$12-\$13) w l lt 3 lw 2 title "He(0)", \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):(\$14/1.0e5) w l lt 0 lw 2 title "T(10^5K)", \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):12 w p lt 1 pt 1 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):13 w p lt 2 pt 2 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):(1-\$12-\$13) w p lt 3 pt 3 notitle
#, \
#     '${f1}' u (${offset}-\$1/1.0e15-0.5):(\$14/1.0e5) w lp lt 4 pt 4 notitle

set output "He_res23.eps"
plot '${f3}' u (${offset}-\$1/1.0e15):12 w l lt 1 lw 2 title "He(1+)", \
     '${f3}' u (${offset}-\$1/1.0e15):13 w l lt 2 lw 2 title "He(2+)", \
     '${f3}' u (${offset}-\$1/1.0e15):(1-\$12-\$13) w l lt 3 lw 2 title "He(0)", \
     '${f3}' u (${offset}-\$1/1.0e15):(\$14/1.0e5) w l lt 0 lw 2 title "T(10^5K)", \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):12 w p lt 1 pt 1 notitle, \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):13 w p lt 2 pt 2 notitle, \
     '${f2}' u (${offset}-\$1/1.0e15-0.25):(1-\$12-\$13) w p lt 3 pt 3 notitle
#, \
#    '${f2}' u (${offset}-\$1/1.0e15-0.25):(\$14/1.0e5) w lp lt 4 pt 4 notitle

set output "He_res14.eps"
plot '${f4}' u (${offset}-\$1/1.0e15+0.1):12 w l lt 1 lw 2 title "He(1+)", \
     '${f4}' u (${offset}-\$1/1.0e15+0.1):13 w l lt 2 lw 2 title "He(2+)", \
     '${f4}' u (${offset}-\$1/1.0e15+0.1):(1-\$12-\$13) w l lt 3 lw 2 title "He(0)", \
     '${f4}' u (${offset}-\$1/1.0e15+0.1):(\$14/1.0e5) w l lt 0 lw 2 title "T(10^5K)", \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):12 w p lt 1 pt 1 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):13 w p lt 2 pt 2 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):(1-\$12-\$13) w p lt 3 pt 3 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):(\$14/1.0e5) w l lt 0 lw 2 lc 4 notitle

set output "He_res13.eps"
plot '${f3}' u (${offset}-\$1/1.0e15):12 w l lt 1 lw 2 title "He(1+)", \
     '${f3}' u (${offset}-\$1/1.0e15):13 w l lt 2 lw 2 title "He(2+)", \
     '${f3}' u (${offset}-\$1/1.0e15):(1-\$12-\$13) w l lt 3 lw 2 title "He(0)", \
     '${f3}' u (${offset}-\$1/1.0e15):(\$14/1.0e5) w l lt 0 lw 2 title "T(10^5K)", \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):12 w p lt 1 pt 1 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):13 w p lt 2 pt 2 notitle, \
     '${f1}' u (${offset}-\$1/1.0e15-0.75):(1-\$12-\$13) w p lt 3 pt 3 notitle
#, \
#     '${f1}' u (${offset}-\$1/1.0e15-0.75):(\$14/1.0e5) w lp lt 4 pt 4 notitle
EOF
gnuplot gnu2.plt
exit

simfile=$1
cat << EOF  > gnu.plt
set size 1.0,0.75
set xrange [2.4e16:3.4e16]
set yrange [-0.1:1.5]
#set yrange [-0.1:1.1]
#set yrange [0.9e-23:1.5e-21]
#set log y
EOF

LIST=`ls ${simfile}*`
for i in $LIST
do
    cat << EOF  >> gnu.plt
#plot '${i}' u 1:2 w l lt -1 
#plot '${i}' u 1:10 w l title "x(e-)", \
# '${i}' u 1:11 w l title "x(H+)", \
# '${i}' u 1:(1-\$11) w l title "x(H0)", \
# '${i}' u 1:(\$14/1.0e5) w l lt -1 title "T(10^5K)"
plot '${i}' u 1:12 w l title "He(1+)", \
 '${i}' u 1:13 w l title "He(2+)", \
 '${i}' u 1:(1-\$12-\$13) w l title "He(0)", \
 '${i}' u 1:(\$14/1.0e5) w l lt -1 title "T(10^5K)"
pause 0.1
EOF
done

gnuplot gnu.plt
