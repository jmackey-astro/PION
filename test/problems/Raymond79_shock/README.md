Shock test based on [Raymond (1979)](https://ui.adsabs.harvard.edu/abs/1979ApJS...39....1R/abstract), model E.

* The left boundary is a reflecting wall, and initial conditions are a flow at 100 km/s to the left.
* The right boundary is in inflow condition.


Commands to run the simulation (example):

```
# Run simulation
../../../build/icgen-ug pf_RShRay79_wB_HHe_ModelE_n0256.txt silo
../../../build/pion-ug RSH1D_n0256_v100_Ray79E_HHe_0000.00000000.silo solver=8

mkdir -p img
# plot images of each snapshot:
python3 plot_1D_state.py ./ RSH1D_n0256_v100_Ray79E_HHe_0000 img

```

If you use gnuplot, you can plot a snapshot with e.g.
```
# Generate text files with the simulation data vs. position:
../../../build/silo2text ./ RSH1D_n0256_v100_Ray79E_HHe_0000. txt_RSH1D_n0256_v100_Ray79E_HHe 1
gnuplot
> # plot Nitrogen ionization states:
> plot "txt_RSH1D_n0256_v100_Ray79E_HHe.00040960.txt" u 1:24 w l title "N+", "" u 1:25 w l title "N2+", "" u 1:26 w l title "N3+", "" u 1:27 w l title "N4+", "" u 1:($13-($24+$25+$26+$27+$28)) w lp title "N0", "" u 1:13 w l lw 2 title "N"
```

