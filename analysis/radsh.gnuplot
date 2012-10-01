set term postscript enhanced colour # eps
#

set output "RSh1d_Raymond1979E_633yr_rpvt.eps"
set size 1.0,1.0
set origin 0.0,0.0
set label 1 "100km/s Magnetised Radiative Shock Properties (t=633yr)" at screen 0.2, screen 0.98
show label 1
unset key 
#
set multiplot
#
set xrange [1.2e16:0]
unset key
set lmargin 15
set rmargin 1
# first plot
set size 0.9, 0.27
set origin 0.0,0.0
set xlabel "x (cm)" 0,0.2
set ylabel "Density" 0.7
set key noautotitles
set log y
plot 'higherP.106328'  using 1:2 w l
#
# SECOND PLOT
set origin 0.0,0.23
set size 0.9, 0.27
set yrange [*:*]
set log y
#set yrange [1.e-12:*]
set format x ""
set xlabel ""
set ylabel "Pressure" 0.7
plot 'higherP.106328'  using 1:3 with l 
#
# THIRD PLOT
set origin 0.0,0.46
unset log y
set yrange [*:*]
set ylabel "Velocity" 0.7
plot 'higherP.106328'  using 1:4 w l 
#
# FOURTH PLOT
set origin 0.0,0.69
set log y
set yrange [*:*]
set ylabel "T (K)" 0.7
plot 'higherP.106328'  using 1:15 w l
#
unset multiplot
#
#
pause -1

set size 1.0,1.0
set origin 0.0,0.0
unset lmargin
unset rmargin
unset log y
unset label 1
set format x "%g"
set output "RSh1d_Raymond1979E_633yr_dens.eps"
set xrange [1.25e16:0.]
set xlabel "Position [cm]"
set ylabel "Density [cm^{-3}]/Temperature [10^2K]"
set title "Density/Temperature through 100km/s Radiative Shock (B_0=1muG)"
set ytics 100
set grid
#set log y
#set yrange [1:1000]
plot "higherP.106328" u 1:($2*$11/1.38/1.67e-24) w l lw 3 title "Electron Density", \
     "higherP.106328" u 1:($2/1.26/1.67e-24) w l lw 3     title "  Nuclei Density", \
     "higherP.106328" u 1:($15/100) w l lw 3              title "       T (10^2K)"
pause -1

set output "RSh1d_Raymond1979E_633yr_H.eps"
set xrange [1.25e16:0.]
set xlabel "Position [cm]"
set ylabel "Concentration/Temperature"
set title "Hydrogen Ionisation through 100km/s Radiative Shock"
unset log y
set ytics 0.1
set grid
#set log y
#set yrange [1.e-2:1.5]
plot "higherP.106328" u 1:(1-$12) w l lw 3 title "H I  ", \
     "higherP.106328" u 1:12 w l lw 3 title "H II ", \
     "higherP.106328" u 1:11 w l lw 3 title "elec.", \
     "higherP.106328" u 1:($15/1.e5) w l lw 3 title "T (10^5K)"
pause -1

set output "RSh1d_Raymond1979E_633yr_He.eps"
set xrange [1.25e16:0.]
set xlabel "Position [cm]"
set ylabel "Concentration/Temperature"
set title "Helium Ionisation through 100km/s Radiative Shock"
unset log y
set ytics 0.1
set grid
#set log y
#set yrange [1.e-2:1.5]
plot "higherP.106328" u 1:(1-$13-$14) w l lw 3 title "He I  ", \
     "higherP.106328" u 1:13 w l lw 3 title "He II ", \
     "higherP.106328" u 1:14 w l lw 3 title "He III", \
     "higherP.106328" u 1:($15/1.e5) w l lw 3 title "T (10^5K)"
pause -1
