#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import sys
sys.path.append("/home/jm/active/projects/silo/test/lib")
import Silo
sys.path.append("/home/jm/code/pypion/Library")
import Plotting_Classes as ppion

for file in ["advection_v020_t30_level00_0000.00000000.silo","advection_v020_t30_level00_0000.00007394.silo", "adv_l1_0000.00003697.silo"]:
  plt.figure()
  plotting=ppion.Plotting2d(file)
  n = plotting.nlevels()
  c = plotting.cycle()
  t = plotting.sim_time().value

  for ilev in range(n):
    #if ilev==0 and n>1:
    #  continue
    if (n>1):
      s1 = "level"+ "%02d" % ilev
      if (file.find(s1) == -1):
        last=ilev-1
        s2 = "level"+ "%02d" % last
        file=file.replace(s2,s1)
        if (file.find(s1) == -1):
          print("error with filename, f=",file)
          quit()

    plotting=ppion.Plotting2d(file)
    #d = plotting.reshaped_parameter2d("VelocityX")
    d = plotting.reshaped_parameter2d("Density")
    level_min = plotting.level_min()
    level_max = plotting.level_max()
    level_ng  = plotting.ngrid()
    dy=((level_max[1]-level_min[1])).value/level_ng[1]
    dx=((level_max[0]-level_min[0])).value/level_ng[0]
    y, x = np.mgrid[slice(level_min[1].value, level_max[1].value, dy),
                    slice(level_min[0].value, level_max[0].value, dx)]
    #print dy, dx
    #print y, x
    if (ilev==0):
      plt.xlim(level_min[0].value, level_max[0].value)
      plt.ylim(level_min[1].value, level_max[1].value)
      plt.xlim(0.5,1.5)
      plt.ylim(0.5,1.5)
    
      #plt.pcolormesh(x, y, d, cmap="viridis",linewidth=0,rasterized=True)
      #clev = MaxNLocator(nbins=7).tick_values(d.min(), d.max())
      clev = [0.99999,1.005, 1.1]
      #print(clev, d.min(), d.max())
      plt.contour(x, y, d, colors="black",levels=clev,linewidths=0.5)

    clev = MaxNLocator(nbins=11).tick_values(0.5,10.5)
    clev = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]
    #print(clev)
    plt.contourf(x, y, d, cmap="viridis",levels=clev,linewidths=0,rasterized=True)

    #extents = [level_min[0].value, level_max[0].value, \
    #           level_min[1].value,level_max[1].value]
    #plt.imshow(d,interpolation="nearest",cmap="magma", \
    #           extent=extents,origin="lower", \
    #           vmin=0.0,vmax=10.0)
    
    plotting.close()
    del plotting

  plt.colorbar()
  plt.xlabel("x")
  plt.ylabel("y")
#plt.grid()
  dmax=str("%.6f" % (d.max()-10.0))
  dmin=str("%.6f" % (d.min()-1.0))
  tm = str("%.5f" % t)
  st = "Density, t=2, max=" +dmax + " min=" + dmin
  st = "Density, t="+tm
  plt.text(0.51,1.43,st,color="white")
#plt.show()
  fn = "advection2D_l"+str(n)+"c"+str(c)+".png"
  plt.savefig(fn,bbox_inches="tight",dpi=300)
quit()
