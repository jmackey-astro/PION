#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/jm/active/projects/silo/test/lib")
import Silo
sys.path.append("/home/jm/code/pypion/Library")
import Plotting_Classes as ppion
import glob

solver="HLL"
solver="Roe"
if solver=="HLL":
  files=["HLL_v020t30_l1n128_0000.00001847.silo", 
             "HLL_v020t30_l1n256_0000.00003694.silo", 
             "HLL_v020t30_l1n512_0000.00000000.silo", 
             "HLL_v020t30_l1n512_0000.00007388.silo", 
             "HLL_v020t30_l2n128_level00_0000.00003700.silo",
             "HLL_v020t30_l2n256_level00_0000.00007394.silo",
             "HLL_v020t30_l2n512_level00_0000.00014782.silo"]
else:
  files=["advection_v020_t30_l1n128_0000.00001847.silo", 
             "advection_v020_t30_l1n256_0000.00003694.silo", 
             "advection_v020_t30_l1n512_0000.00007388.silo", 
             "advection_v020_t30_l2n128_level00_0000.00003700.silo",
             "advection_v020_t30_l2n256_level00_0000.00007394.silo",
             "advection_v020_t30_l2n512_level00_0000.00014782.silo"]

f0=sorted(glob.glob("april/advection_v020_t30_l2n128_level00_0000.*.silo"))
f1=sorted(glob.glob("april/advection_v020_t30_l2n128_level01_0000.*.silo"))
f0=sorted(glob.glob("advection_v020_t30_l2n*_level00_0000.*.silo"))
f1=sorted(glob.glob("advection_v020_t30_l2n*_level01_0000.*.silo"))
#f0=sorted(glob.glob("advection_v020_t30_l1n128_0000.*.silo"))
#f0=sorted(glob.glob("advection_v020_t30_l1n128.*.silo"))
print(f0)
for i in range(0,len(f0)):
  file = [f0[i],f1[i]]
  #file = [f0[i]]
  plt.figure()
  plotting=ppion.Plotting2d(file)
  n = plotting.nlevels()
  c = plotting.cycle()
  cy = str(c).zfill(8)
  t = plotting.sim_time().value
  level_ng  = plotting.ngrid()
  D = plotting.get_2Darray("Density")
  d = D['data']
  t = D['sim_time']
  print("time=",t)
  xmax = D['max_extents']
  xmin = D['min_extents']

  for ilev in range(n):
    level_min = xmin[ilev]
    level_max = xmax[ilev]
    dy=((level_max[1]-level_min[1]))/level_ng[1]
    dx=((level_max[0]-level_min[0]))/level_ng[0]
    y, x = np.mgrid[slice(level_min[1]+0.5*dy, level_max[1]+0.5*dy, dy),
                    slice(level_min[0]+0.5*dx, level_max[0]+0.5*dx, dx)]

    if (ilev==0):
      plt.xlim(level_min[0], level_max[0])
      plt.ylim(level_min[1], level_max[1])
      plt.xlim(0.4,1.6)
      plt.ylim(0.4,1.6)
    
      #plt.pcolormesh(x, y, d, cmap="viridis",linewidth=0,rasterized=True)
      #clev = MaxNLocator(nbins=7).tick_values(d.min(), d.max())
      clev = [0.99999,1.005, 1.1]
      #print(clev, d.min(), d.max())
      plt.contour(x, y, d[ilev], colors="white",levels=clev,linewidths=0.5)
    #else:
    #clev = MaxNLocator(nbins=11).tick_values(0.5,10.5)
    clev = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]
    #print(clev)
    plt.contourf(x, y, d[ilev], cmap="viridis",levels=clev)

    #extents = [level_min[0].value, level_max[0].value, \
    #           level_min[1].value,level_max[1].value]
    #plt.imshow(d,interpolation="nearest",cmap="magma", \
    #           extent=extents,origin="lower", \
    #           vmin=0.0,vmax=10.0)
    
  plt.colorbar()
  plt.xlabel("x",fontsize=16, fontweight='bold')
  plt.ylabel("y",fontsize=16, fontweight='bold')
  #plt.grid()
  dmax=str("%.6f" % (d[0].max()))
  dmin=str("%.6f" % (d[0].min()))
  tm = str("%.5f" % t)
  st = "Density, max=" +dmax + " min=" + dmin
  #st = "Density, t="+tm
  plt.text(0.41,1.53,st,color="white")
  #plt.show()

  fn = file[0][0:25] + "c"+cy+".png"
  #fn = file[0][6:31] + "c"+cy+".png"
  plt.savefig(fn,bbox_inches="tight",dpi=300)

quit()
