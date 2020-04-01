#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/jm/active/projects/silo/test/lib")
import Silo
sys.path.append("/home/jm/code/pypion/Library")
import Plotting_Classes as ppion

#plt.rc('font',**{'family':'sans-serif','sans-serif':['Bitstream Vera Sans']})
#plt.rc('text', usetex=True)
#plt.rc('font',**{'size': 14})
#plt.rc('lines', linewidth=2)
#plt.rcParams['mathtext.fontset'] = 'stixsans'
#plt.rcParams['font.family'] = 'Bitstream Vera Sans'
#plt.rcParams['mathtext.rm'] = 'stixsans'
#plt.rcParams['mathtext.it'] = 'stixsans:italic'
#plt.rcParams['mathtext.bf'] = 'stixsans:bold'
plt.rcParams["font.weight"] = "normal"
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_family('sans-serif')
font.set_name('stixsans')
font.set_style('italic')
font.set_weight('light')

#solver="HLL"
solver="HLLD"
#solver="Roe"
#solver="RCV"
solver="NG"
if solver=="HLL":
  files=["HLL_n100_0000.00000000.silo","HLL_n100_0000.00000412.silo","HLL_n100_0000.00000825.silo"]
elif solver=="HLLD":
  files=["HLLD_n100_0000.00000000.silo","HLLD_n100_0000.00000825.silo","HLLD_n200_0000.00000000.silo","HLLD_n200_0000.00001646.silo"]
elif solver=="RCV":
  files=["RCV_n100_0000.00000000.silo","RCV_n100_0000.00000412.silo","RCV_n100_0000.00000825.silo","RCV_n200_0000.00000000.silo","RCV_n200_0000.00000823.silo","RCV_n200_0000.00001646.silo","FieldLoop200l2_level00.00000000.silo"]
elif solver=="NG":
  files=["FieldLoop200l2_HLLD_level00_0000.00000000.silo","FieldLoop200l2_HLLD_level00_0000.00003316.silo"]
else:
  files=["FieldLoop100_0000.00000000.silo",  "FieldLoop100_0000.00000412.silo",  "FieldLoop100_0000.00000825.silo"]
  #files=["FieldLoop032_0000.00000000.silo", "FieldLoop032_0000.00000132.silo", "FieldLoop032_0000.00000264.silo"]

for file in files:
  fig, (ax1,ax2) = plt.subplots(2,1)
  
  dataio=ppion.Plotting2d(file)
  n = dataio.nlevels()
  c = dataio.cycle()
  t = dataio.sim_time().value
  level_ng  = dataio.ngrid()

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

    dataio=ppion.Plotting2d(file)
    #vx = dataio.reshaped_parameter2d("VelocityX")
    #d = dataio.reshaped_parameter2d("Density")
    bx = dataio.reshaped_parameter2d("MagneticFieldX")
    by = dataio.reshaped_parameter2d("MagneticFieldY")
    bz = dataio.reshaped_parameter2d("MagneticFieldZ")
    cB = dataio.reshaped_parameter2d("CurlB")
    B = np.sqrt(4.0*np.pi*(bx*bx + by*by + bz*bz))
    Pmag=B**2/2.0
    
    level_min = dataio.level_min()
    level_max = dataio.level_max()
    dy=((level_max[1]-level_min[1])).value/level_ng[1]
    dx=((level_max[0]-level_min[0])).value/level_ng[0]
    y, x = np.mgrid[slice(level_min[1].value+0.5*dy, level_max[1].value+0.5*dy, dy),
                    slice(level_min[0].value+0.5*dx, level_max[0].value+0.5*dx, dx)]
    if (ilev==0):
      ax1.set_xlim(level_min[0].value, level_max[0].value)
      ax1.set_ylim(level_min[1].value, level_max[1].value)
      ax1.set_xlim(-1.0,1.0)
      ax1.set_ylim(-0.5,0.5)
      ax1.set_ylabel("y",fontsize=14, fontweight='bold')
      
      ax2.set_xlim(level_min[0].value, level_max[0].value)
      ax2.set_ylim(level_min[1].value, level_max[1].value)
      ax2.set_xlim(-1.0,1.0)
      ax2.set_ylim(-0.5,0.5)
      ax2.set_xlabel("x",fontsize=14, fontweight='bold')
      ax2.set_ylabel("y",fontsize=14, fontweight='bold')
    if (ilev>0) || (n==1):
      #plt.pcolormesh(x, y, d, cmap="viridis",linewidth=0,rasterized=True)
      clev = MaxNLocator(nbins=7).tick_values(B.min(), B.max())
      #plt.contour(x, y, B, colors="white",levels=clev,linewidths=0.5,aspect=1)
      clev = [0.0,0.002,0.004,0.006,0.008,0.010,0.012,0.014]
      print clev
      ct1 = ax1.contour(x, y, B, colors="white",levels=clev,linewidths=0.5)
      ct1.monochrome = True

      clev = MaxNLocator(nbins=7).tick_values(cB.min(), cB.max())
      #clev = [-0.016 -0.008  0.     0.008  0.016  0.024  0.032  0.04 ]
      print clev
      ct2 = ax2.contour(x, y, cB, colors="white",levels=clev,linewidths=0.5)
      ct2.monochrome = True


    #clev = MaxNLocator(nbins=11).tick_values(d.min(), d.max())
    #print(clev)
    #plt.contourf(x, y, B, cmap="viridis",levels=clev)
    #plt.contourf(x, y, cb, cmap="magma",levels=clev)

    extents = [level_min[0].value, level_max[0].value, \
               level_min[1].value,level_max[1].value]
    #plt.imshow(B,interpolation="nearest",cmap="viridis", \
    #           extent=extents,origin="lower")
    vn=0.0; vp=0.013
    im1 = ax1.imshow(B,interpolation="nearest",cmap="viridis", \
               extent=extents,origin="lower",vmin=vn,vmax=vp)
    if c==0:
      vn=-0.13; vp=0.275
    else:
      vn=-0.015; vp=0.037
    im2 = ax2.imshow(cB,interpolation="nearest",cmap="viridis", \
               extent=extents,origin="lower",vmin=vn,vmax=vp)

    dataio.close()
    del dataio

  fig.colorbar(im1, ax=ax1)
  fig.colorbar(im2, ax=ax2)

  dmax=str("%.4f" % (B.max()))
  dmin=str("%.4f" % (B.min()))
  tm = str("%.5f" % t)
  #st = "$P_\mathrm{mag}$, max=" +dmax + " min=" + dmin
  st = "$\left| B \\right|$, max=" +dmax + " min=" + dmin
  #st = "Density, t="+tm
  ax1.text(-0.75,0.42,st,color="white",fontsize=10)
  if c==0:
    ax1.text(-0.95,0.42,"(a)",color="white",fontsize=10,weight='bold')
  else:
    ax1.text(-0.95,0.42,"(c)",color="white",fontsize=10,weight='bold')

  dmax=str("%.4f" % (cB.max()))
  dmin=str("%.4f" % (cB.min()))
  tm = str("%.5f" % t)
  st = "$\\nabla\\times B$, max=" +dmax + " min=" + dmin
  #st = "Density, t="+tm
  ax2.text(-0.75,0.42,st,color="black",fontsize=10)
  if c==0:
    ax2.text(-0.95,0.42,"(b)",color="black",fontsize=10,weight='bold')
  else:
    ax2.text(-0.95,0.42,"(d)",color="black",fontsize=10,weight='bold')

  for ax in fig.get_axes():
    ax.label_outer()

  #plt.show()
  fn = file[0:12] + "_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".png"
  plt.savefig(fn,bbox_inches="tight",dpi=300)

#  if solver=="HLL":
#    fn = "HLL2D_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".png"
#    plt.savefig(fn,bbox_inches="tight",dpi=300)
#    fn = "HLL2D_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".pdf"
#    plt.savefig(fn,bbox_inches="tight",dpi=300)
#  elif solver=="RCV":
#    fn = "RCV2D_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".png"
#    plt.savefig(fn,bbox_inches="tight",dpi=300)
#    fn = "RCV2D_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".pdf"
#    plt.savefig(fn,bbox_inches="tight",dpi=300)
#  else:
#    fn = "HLLD2D_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".png"
#    plt.savefig(fn,bbox_inches="tight",dpi=300)
#    fn = "HLLD2D_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".pdf"
#    plt.savefig(fn,bbox_inches="tight",dpi=300)
quit()
