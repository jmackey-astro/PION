
set terminal wxt enhanced dashed size 800,600 font "Arial,20"

set lmargin 8.0
set rmargin 1.0
set tmargin 1.0
set bmargin 3.0
set border linewidth 2.0

set style line 11 lt 1 lc rgb "red" lw 3 pt 1 ps 1.5
set style line 12 lt 2 lc rgb "red" lw 3 pt 2 ps 1.5
set style line 13 lt 3 lc rgb "red" lw 3 pt 3 ps 1.5
set style line 14 lt 4 lc rgb "red" lw 3 pt 4 ps 1.5
set style line 15 lt 5 lc rgb "red" lw 3 pt 5 ps 1.5
set style line 16 lt 6 lc rgb "red" lw 3 pt 6 ps 1.5

set style line 21 lt 1 lc rgb "gray20" lw 3 pt 1 ps 1.5
set style line 22 lt 2 lc rgb "gray20" lw 3 pt 2 ps 1.5
set style line 23 lt 3 lc rgb "gray20" lw 3 pt 3 ps 1.5
set style line 24 lt 4 lc rgb "gray20" lw 3 pt 4 ps 1.5
set style line 25 lt 5 lc rgb "gray20" lw 3 pt 5 ps 1.5
set style line 26 lt 6 lc rgb "gray20" lw 3 pt 6 ps 1.5

set style line 31 lt 1 lc rgb "blue" lw 3 pt 1 ps 1.5
set style line 32 lt 2 lc rgb "blue" lw 3 pt 2 ps 1.5
set style line 33 lt 3 lc rgb "blue" lw 3 pt 3 ps 1.5
set style line 34 lt 4 lc rgb "blue" lw 3 pt 4 ps 1.5
set style line 35 lt 5 lc rgb "blue" lw 3 pt 5 ps 1.5
set style line 36 lt 6 lc rgb "blue" lw 3 pt 6 ps 1.5

set style line 41 lt 1 lc rgb "web-green" lw 3 pt 1 ps 1.5
set style line 42 lt 2 lc rgb "web-green" lw 3 pt 2 ps 1.5
set style line 43 lt 3 lc rgb "web-green" lw 3 pt 3 ps 1.5
set style line 44 lt 4 lc rgb "web-green" lw 3 pt 4 ps 1.5
set style line 45 lt 5 lc rgb "web-green" lw 3 pt 5 ps 1.5
set style line 46 lt 6 lc rgb "web-green" lw 3 pt 6 ps 1.5


set style line 55 lt 2 lc rgb "gray20" lw 2 pt 2 ps 1.5



h_cr(xe)       = 5.0e-28*(1.0-xe)
h_fuv(xe,T,n)  = 1.083e-25*(1.0-xe)/(1.0+9.77e-3*(sqrt(T)/(n*xe))**0.73)

c_fbd(xe,T,n)  = 1.20e-22*exp(-33610.0/T -(2180.0*2180.0/T/T)) *xe*xe*n*exp(-T*T/5.0e10)
# Added in density dependence of CII line for collisions with H (for consistency with CII-e rate).
#c_CIIH(xe,T,n) = 3.15e-27*exp(-92.0/T)*n*(1.0-xe)*exp(-n/1.0e4/(1.0 + 0.05*n*(T/2000.0)**(-0.37)))
c_CIIH(xe,T,n) = 3.15e-27*exp(-92.0/T)*n*(1.0-xe)*exp(-n/1.0e4)
c_OIH(xe,T,n)  = 3.96e-28*exp(0.4*log(T)-228.0/T)*n*(1.0-xe)
#c_CIIe(xe,T,n) = 1.4e-23*exp(-0.5*log(T)-92.0/T)*n*xe*exp(-n/1.0e4/(1.0 + 0.05*n*(T/2000.0)**(-0.37)))
c_CIIe(xe,T,n) = 1.4e-23*exp(-0.5*log(T)-92.0/T)*n*xe*exp(-n/1.0e4)
c_PAH(xe,T,n)  = 3.02e-30*exp(0.94*log(T) +0.74*T**(-0.068)*log(3.4*sqrt(T)/(n*xe)))*n*xe
c_TOT(xe,T,n) = c_fbd(xe,T,n) + c_CIIH(xe,T,n) + c_OIH(xe,T,n) + c_CIIe(xe,T,n) + c_PAH(xe,T,n)

# Wolfire et al. 2003, approximate rate for cooling in 100-1000K range.
c_WOLF(xe,T,n) = 5.4e-27*exp(-150.0/T +0.2*log(T/100.0))*n*(1.0-xe)*exp(-T*T/5.0e10)

# Koyama and Inutsuka (2002) cooling rate.
c_KI(T,n) = n* (2.0e-19*exp(-1.184e5/(T+1.0e3)) +2.8e-28*sqrt(T)*exp(-92.0/T))

# PDR cooling rate from Henney et al. 2009.
c_MOL(xe,T,n) = 3.981e-27* n**0.6* sqrt(T)*exp(-(70.0 +220.0*(n/1.0e6)**0.2)/T)

# MacDonald & Bailey 1981 cooling rate at T<10^4K (used by Villaver) (Approximate)
#c_MB(T,n) = 4e-24*n*(T*1.0e-4)**2.25 *exp(-T*T/1.0e9)

set xrange [1e1:1e5]
set log x
set grid x
set xlabel "Temperature (K)" offset 0,0.5
set xtics 1.0,10 format  "10^{%T}" offset 0.25,0.2

set yrange [1e-30:1e-21]
set log y
set grid y
set ytics 1.0e-30,10 format  "10^{%T}" offset 0.25,0.2
set ylabel "Rate (erg.cm^3/s)"

set key bottom right

Ke=0.1
NN=1.0
set label 1 "n(H)=1 cm^{-3}, x_e = 0.1" at screen 0.2,0.9
plot  c_fbd(Ke,x,NN)/NN w l ls 11 title "Forbidden", \
      c_CIIH(Ke,x,NN)/NN w l ls 21 title "CII-H coll.", \
      c_OIH(Ke,x,NN)/NN w l ls 22 title "OI-H coll.", \
      c_CIIe(Ke,x,NN)/NN w l ls 23 title "CII-e coll.", \
      c_PAH(Ke,x,NN)/NN w l ls 31 title "PAH coll.", \
      c_TOT(Ke,x,NN)/NN w l ls 41 title "Total", \
      c_KI(x,NN)/NN w l ls 43 title "KI02", \
      c_MOL(Ke,x,NN)/NN w l lt -1 title "Molecular"
#, \
#    c_MB(x,NN)/NN w l lt -1 title "MB81"
#      c_WOLF(Ke,x,NN)/NN w l ls 42 title "Wolfire C4", \

pause -1

Ke=0.1
NN=10.0
set label 1 "n(H)=10.0 cm^{-3}, x_e = 0.1" at screen 0.2,0.9
plot  c_fbd(Ke,x,NN)/NN w l ls 11 title "Forbidden", \
      c_CIIH(Ke,x,NN)/NN w l ls 21 title "CII-H coll.", \
      c_OIH(Ke,x,NN)/NN w l ls 22 title "OI-H coll.", \
      c_CIIe(Ke,x,NN)/NN w l ls 23 title "CII-e coll.", \
      c_PAH(Ke,x,NN)/NN w l ls 31 title "PAH coll.", \
      c_TOT(Ke,x,NN)/NN w l ls 41 title "Total", \
      c_KI(x,NN)/NN w l ls 43 title "KI02", \
      c_MOL(Ke,x,NN)/NN w l lt -1 title "Molecular"
#, \
#      c_MB(x,NN)/NN w l lt -1 title "MB81"
pause -1

Ke=0.1
NN=100.0
set label 1 "n(H)=100.0 cm^{-3}, x_e = 0.1" at screen 0.2,0.9
plot  c_fbd(Ke,x,NN)/NN w l ls 11 title "Forbidden", \
      c_CIIH(Ke,x,NN)/NN w l ls 21 title "CII-H coll.", \
      c_OIH(Ke,x,NN)/NN w l ls 22 title "OI-H coll.", \
      c_CIIe(Ke,x,NN)/NN w l ls 23 title "CII-e coll.", \
      c_PAH(Ke,x,NN)/NN w l ls 31 title "PAH coll.", \
      c_TOT(Ke,x,NN)/NN w l ls 41 title "Total", \
      c_KI(x,NN)/NN w l ls 43 title "KI02", \
      c_MOL(Ke,x,NN)/NN w l lt -1 title "Molecular"
#, \
#      c_MB(x,NN)/NN w l lt -1 title "MB81"
pause -1

Ke=0.01
NN=1.0
set label 1 "n(H)=1 cm^{-3}, x_e = 0.01" at screen 0.2,0.9
plot  c_fbd(Ke,x,NN)/NN w l ls 11 title "Forbidden", \
      c_CIIH(Ke,x,NN)/NN w l ls 21 title "CII-H coll.", \
      c_OIH(Ke,x,NN)/NN w l ls 22 title "OI-H coll.", \
      c_CIIe(Ke,x,NN)/NN w l ls 23 title "CII-e coll.", \
      c_PAH(Ke,x,NN)/NN w l ls 31 title "PAH coll.", \
      c_TOT(Ke,x,NN)/NN w l ls 41 title "Total", \
      c_KI(x,NN)/NN w l ls 43 title "KI02", \
      c_MOL(Ke,x,NN)/NN w l lt -1 title "Molecular"
#, \
#      c_MB(x,NN)/NN w l lt -1 title "MB81"
pause -1

Ke=0.01
NN=10.0
set label 1 "n(H)=10 cm^{-3}, x_e = 0.01" at screen 0.2,0.9
plot  c_fbd(Ke,x,NN)/NN w l ls 11 title "Forbidden", \
      c_CIIH(Ke,x,NN)/NN w l ls 21 title "CII-H coll.", \
      c_OIH(Ke,x,NN)/NN w l ls 22 title "OI-H coll.", \
      c_CIIe(Ke,x,NN)/NN w l ls 23 title "CII-e coll.", \
      c_PAH(Ke,x,NN)/NN w l ls 31 title "PAH coll.", \
      c_TOT(Ke,x,NN)/NN w l ls 41 title "Total", \
      c_KI(x,NN)/NN w l ls 43 title "KI02", \
      c_MOL(Ke,x,NN)/NN w l lt -1 title "Molecular"
#, \
#      c_MB(x,NN)/NN w l lt -1 title "MB81"
pause -1

Ke=0.01
NN=100.0
set label 1 "n(H)=100.0 cm^{-3}, x_e = 0.01" at screen 0.2,0.9
plot  c_fbd(Ke,x,NN)/NN w l ls 11 title "Forbidden", \
      c_CIIH(Ke,x,NN)/NN w l ls 21 title "CII-H coll.", \
      c_OIH(Ke,x,NN)/NN w l ls 22 title "OI-H coll.", \
      c_CIIe(Ke,x,NN)/NN w l ls 23 title "CII-e coll.", \
      c_PAH(Ke,x,NN)/NN w l ls 31 title "PAH coll.", \
      c_TOT(Ke,x,NN)/NN w l ls 41 title "Total", \
      c_KI(x,NN)/NN w l ls 43 title "KI02", \
      c_MOL(Ke,x,NN)/NN w l lt -1 title "Molecular"
#, \
#      c_MB(x,NN)/NN w l lt -1 title "MB81"
pause -1


Ke=0.001
NN=1.0
set label 1 "n(H)=1 cm^{-3}, x_e = 0.001" at screen 0.2,0.9
plot  c_fbd(Ke,x,NN)/NN w l ls 11 title "Forbidden", \
      c_CIIH(Ke,x,NN)/NN w l ls 21 title "CII-H coll.", \
      c_OIH(Ke,x,NN)/NN w l ls 22 title "OI-H coll.", \
      c_CIIe(Ke,x,NN)/NN w l ls 23 title "CII-e coll.", \
      c_PAH(Ke,x,NN)/NN w l ls 31 title "PAH coll.", \
      c_TOT(Ke,x,NN)/NN w l ls 41 title "Total", \
      c_KI(x,NN)/NN w l ls 43 title "KI02", \
      c_MOL(Ke,x,NN)/NN w l lt -1 title "Molecular"
#, \
#      c_MB(x,NN)/NN w l lt -1 title "MB81"
pause -1

Ke=0.001
NN=10.0
set label 1 "n(H)=10 cm^{-3}, x_e = 0.001" at screen 0.2,0.9
plot  c_fbd(Ke,x,NN)/NN w l ls 11 title "Forbidden", \
      c_CIIH(Ke,x,NN)/NN w l ls 21 title "CII-H coll.", \
      c_OIH(Ke,x,NN)/NN w l ls 22 title "OI-H coll.", \
      c_CIIe(Ke,x,NN)/NN w l ls 23 title "CII-e coll.", \
      c_PAH(Ke,x,NN)/NN w l ls 31 title "PAH coll.", \
      c_TOT(Ke,x,NN)/NN w l ls 41 title "Total", \
      c_KI(x,NN)/NN w l ls 43 title "KI02", \
      c_MOL(Ke,x,NN)/NN w l lt -1 title "Molecular"
#, \
#      c_MB(x,NN)/NN w l lt -1 title "MB81"
pause -1

Ke=0.001
NN=100.0
set label 1 "n(H)=100.0 cm^{-3}, x_e = 0.001" at screen 0.2,0.9
plot  c_fbd(Ke,x,NN)/NN w l ls 11 title "Forbidden", \
      c_CIIH(Ke,x,NN)/NN w l ls 21 title "CII-H coll.", \
      c_OIH(Ke,x,NN)/NN w l ls 22 title "OI-H coll.", \
      c_CIIe(Ke,x,NN)/NN w l ls 23 title "CII-e coll.", \
      c_PAH(Ke,x,NN)/NN w l ls 31 title "PAH coll.", \
      c_TOT(Ke,x,NN)/NN w l ls 41 title "Total", \
      c_KI(x,NN)/NN w l ls 43 title "KI02", \
      c_MOL(Ke,x,NN)/NN w l lt -1 title "Molecular"
#, \
#      c_MB(x,NN)/NN w l lt -1 title "MB81"
pause -1

quit

