
set xrange [100:1e6]
set yrange [1e-30:1e-20]
set log xy
FF_H(x)=1.4e-27*sqrt(x)*ni*ne
FF_He(x)=1.68e-27*0.1*sqrt(x)*ni*ne
FBDN(x)=1.42e-22*exp(-33610.0/x -(2180.0*2180.0/x/x)) *ni*ne*exp(-x*x/5.0e10)
CIIH(x)= 3.15e-27*exp(-92.0/x)*nn*nH
OIH(x)=3.96e-28*exp(0.4*log(x)-228.0/x)*nn*nH
CIIe(x)= 1.4e-23*exp(-0.5*log(x)-92.0/x)*ne*nH
PAH(x)=3.02e-30*exp(0.94*log(x) +0.74*x**(-0.068)*log(3.4*sqrt(x)/ne))*ne*nH
CEXH(x) = 7.5e-19*exp(-118348.0/x)*ne*nn/(1.0+sqrt(1.0e-5*x))
TOTAL(x)=FF_H(x)+ FF_He(x)+ FBDN(x)+ CIIH(x)+ OIH(x)+ CIIe(x)+ PAH(x) +CEXH(x)

nH=1.0; ni=0.999*nH; ne=1.1*ni; nn=nH-ni;
plot FF_H(x), FF_He(x), FBDN(x), CIIH(x), OIH(x), CIIe(x), PAH(x), CEXH(x), TOTAL(x) w l lt -1
pause -1
nH=1.0; ni=0.1*nH; ne=1.1*ni; nn=nH-ni;
plot FF_H(x), FF_He(x), FBDN(x), CIIH(x), OIH(x), CIIe(x), PAH(x), CEXH(x), TOTAL(x) w l lt -1
pause -1
quit


