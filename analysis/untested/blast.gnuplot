set term postscript enhanced color
set output "blast_profile2.eps"
set yrange [1.e-2:5]
set xrange [0.7:30.0]
#set title "Grid Point Data for Blast Wave Test Problem (3d, 128^3)"
set xlabel "Radius (parsecs)"
set ylabel "Pressure, Density, Velocity (scaled)"
#set log y
#set log x
set key left
plot "bw.txt" u 1:2 w p ps 0 title "Density (red)", \
     "bw.txt" u 1:($3*8.e9) w p ps 0 title "Pressure (green)", \
     "bw.txt" u 1:($4/4.e6) w p ps 0 title "Total Velocity (blue)"
pause -1



