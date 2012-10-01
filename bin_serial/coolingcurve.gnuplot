#set term postscript enhanced colour
#set output "cooling_SD93cie_fbdnline.eps"
set ylabel "Cooling Rate (erg/cm^3/s)"
set xlabel "Temperature (K)"
set log y
set log x
set yrange [1.e-28:1.e-10]
#set yrange [1.e-28:1.e-18]
set xrange [1:1.0e7]
set grid
set key top left
#set title "Cooling Rates at n=1.0e4 per cc"
set title "Cooling Rates at n=1.0 per cc, c15"
plot "coolingcurve_c15n4.txt" u 1:2 w lp title "Neutral Gas x=1e-6", \
     "" u 1:3 w lp title "Ion Fraction 0.001", \
     "" u 1:4 w lp title "Ion Fraction 0.1", \
     "" u 1:5 w lp title "Ion Fraction 0.5", \
     "" u 1:6 w lp title "Ion Fraction 0.99", \
     "" u 1:7 w lp title "Ionised Gas  0.999999"
#pause -1
set title "Cooling Rates at n=1.0 per cc, c16"
plot "coolingcurve_c16n4.txt" u 1:2 w lp title "Neutral Gas x=1e-6", \
     "" u 1:3 w lp title "Ion Fraction 0.001", \
     "" u 1:4 w lp title "Ion Fraction 0.1", \
     "" u 1:5 w lp title "Ion Fraction 0.5", \
     "" u 1:6 w lp title "Ion Fraction 0.99", \
     "" u 1:7 w lp title "Ionised Gas  0.999999"
#pause -1
set yrange [1.e-28:1.e-18]
set title "Cooling Rates at n=1.0 per cc, current cooling function"
plot "coolingcurve.txt" u 1:2 w lp title "Neutral Gas x=1e-6", \
     "" u 1:3 w lp title "Ion Fraction 0.001", \
     "" u 1:4 w lp title "Ion Fraction 0.1", \
     "" u 1:5 w lp title "Ion Fraction 0.5", \
     "" u 1:6 w lp title "Ion Fraction 0.99", \
     "" u 1:7 w lp title "Ionised Gas  0.999999"
pause -1
set yrange [1.e-21:1.e-9]
set xrange [1.0:1.e6]
set title "Cooling Rates at n=1.0e6 per cc"
plot "coolingcurve.txt" u 1:8 w lp title "Neutral Gas x=1e-6", \
     "" u 1:9 w lp title "Ion Fraction 0.001", \
     "" u 1:10 w lp title "Ion Fraction 0.1", \
     "" u 1:11 w lp title "Ion Fraction 0.5", \
     "" u 1:12 w lp title "Ion Fraction 0.99", \
     "" u 1:13 w lp title "Ionised Gas  0.999999"
pause -1
quit
