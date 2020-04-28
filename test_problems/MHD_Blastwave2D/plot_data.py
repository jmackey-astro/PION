#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

from matplotlib.ticker import MaxNLocator
from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/jm/active/projects/silo/test/lib")
import Silo
sys.path.append("/home/jm/code/pypion/Library")
import Plotting_Classes as ppion
import glob

plt.rcParams["font.weight"] = "normal"
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_family('sans-serif')
font.set_name('stixsans')
font.set_style('italic')
font.set_weight('light')
font.set_size(8)

solver="HLLD_NG"
#solver="HLLD_UG"
f0=[]
f1=[]
if solver=="HLLD_NG":
  f0=sorted(glob.glob("NG_B*_n256_level00_0000.*.silo"))
  f1=sorted(glob.glob("NG_B*_n256_level01_0000.*.silo"))
elif solver=="HLLD_UG":
  f0=sorted(glob.glob("UG_B*_n256_0000.*.silo"))
else:
  f0=[""]

for i in range(0,len(f0)):
  if len(f1)==len(f0):
    file = [f0[i],f1[i]]
  else:
    file = [f0[i]]
  print(file)
  fig, (ax1,ax2,ax3) = plt.subplots(1,3)
  
  dataio=ppion.Plotting2d(file)
  n = dataio.nlevels()
  c = dataio.cycle()
  level_ng  = dataio.ngrid()

  dataio=ppion.Plotting2d(file)
  D = dataio.get_2Darray("Density")
  d = D['data']
  t = D['sim_time']
  print("time=",t)
  #d = dataio.get_2Darray("Pressure")['data']
  vx = dataio.get_2Darray("VelocityX")['data']
  vy = dataio.get_2Darray("VelocityY")['data']
  vz = dataio.get_2Darray("VelocityZ")['data']
  bx = dataio.get_2Darray("MagneticFieldX")['data']
  by = dataio.get_2Darray("MagneticFieldY")['data']
  bz = dataio.get_2Darray("MagneticFieldZ")['data']
  
  xmax = D['max_extents']
  xmin = D['min_extents']

  for ilev in range(n):
    print("ilev=",ilev)
    level_min = xmin[ilev]
    level_max = xmax[ilev]
    dy=((level_max[1]-level_min[1]))/level_ng[1]
    dx=((level_max[0]-level_min[0]))/level_ng[0]
    y, x = np.mgrid[slice(level_min[1]+0.5*dy, level_max[1]+0.5*dy, dy),
                    slice(level_min[0]+0.5*dx, level_max[0]+0.5*dx, dx)]
    BX=np.array(bx[ilev])
    B = np.sqrt((bx[ilev]*bx[ilev] + by[ilev]*by[ilev] + bz[ilev]*bz[ilev]))
    V2 = 0.5*np.sqrt(vx[ilev]*vx[ilev] + vy[ilev]*vy[ilev] + vz[ilev]*vz[ilev])
    d = np.array(d)
    B = np.array(B)
    Pmag=B**2/(8.0*np.pi)


    if (ilev==0):
      ax1.set_xlim(level_min[0], level_max[0])
      ax1.set_ylim(level_min[1], level_max[1])
      #ax1.set_xlim(-1.0,1.0)
      #ax1.set_ylim(-0.5,0.5)
      ax1.set_xlabel("x",fontsize=10, fontweight='bold')
      ax1.set_ylabel("y",fontsize=10, fontweight='bold')
      ax1.xaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
      ax1.tick_params(labelsize=8)
      
      ax2.set_xlim(level_min[0], level_max[0])
      ax2.set_ylim(level_min[1], level_max[1])
      #ax2.set_xlim(-1.0,1.0)
      #ax2.set_ylim(-0.5,0.5)
      ax2.set_xlabel("x",fontsize=10, fontweight='bold')
      #ax2.set_ylabel("y",fontsize=12, fontweight='bold')
      ax2.yaxis.set_ticklabels([])
      ax2.xaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
      ax2.tick_params(labelsize=8)
      
      ax3.set_xlim(level_min[0], level_max[0])
      ax3.set_ylim(level_min[1], level_max[1])
      #ax3.set_xlim(-1.0,1.0)
      #ax3.set_ylim(-0.5,0.5)
      ax3.set_xlabel("x",fontsize=10, fontweight='bold')
      #ax3.set_ylabel("y",fontsize=12, fontweight='bold')
      ax3.yaxis.set_ticklabels([])
      ax3.xaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
      ax3.tick_params(labelsize=8)
      
    #if (ilev>0) or (n==1):
      #plt.pcolormesh(x, y, d, cmap="viridis",linewidth=0,rasterized=True)
      clev = MaxNLocator(nbins=21).tick_values(d[ilev].min(), d[ilev].max())
      #plt.contour(x, y, B, colors="white",levels=clev,linewidths=0.5,aspect=1)
      #clev = [0.0,0.002,0.004,0.006,0.008,0.010,0.012,0.014]
      ct1 = ax1.contour(x, y, d[ilev], colors="white",levels=clev,linewidths=0.25)
      ct1.monochrome = True

      clev = MaxNLocator(nbins=7).tick_values(V2.min(), V2.max())
      #clev = [-0.016 -0.008  0.     0.008  0.016  0.024  0.032  0.04 ]
      ct2 = ax2.contour(x, y, V2, colors="white",levels=clev,linewidths=0.5)
      ct2.monochrome = True

      clev = MaxNLocator(nbins=7).tick_values(Pmag.min(), Pmag.max())
      #clev = [-0.016 -0.008  0.     0.008  0.016  0.024  0.032  0.04 ]
      ct3 = ax3.contour(x, y, Pmag, colors="white",levels=clev,linewidths=0.5)
      ct3.monochrome = True

    extents = [level_min[0], level_max[0], \
               level_min[1],level_max[1]]

    #clev = MaxNLocator(nbins=11).tick_values(d.min(),d.max())
    #print(clev)
    #im1 =ax1.contourf(x, y, d, cmap="viridis",levels=clev)
    if ilev==0:
      dn=0.0; dp=d[ilev].max()
    im1 = ax1.imshow(d[ilev],interpolation="nearest",cmap="viridis", \
               extent=extents,origin="lower",vmin=dn,vmax=dp)
    if (ilev==0):
      #fig.colorbar(im1, ax=ax1,orientation="horizontal", pad=-0.5)
      cbaxes = fig.add_axes([0.12, 0.8, 0.23, 0.02]) 
      cb=fig.colorbar(im1, ax=ax1,orientation="horizontal", cax=cbaxes,pad=0.0)
      cb.ax.tick_params(labelsize=8)
      tick_locator = ticker.MaxNLocator(nbins=4)
      cb.locator = tick_locator
      cb.update_ticks()

    if ilev==0:
      v2n=0.0; v2p=V2.max()
    im2 = ax2.imshow(V2,interpolation="nearest",cmap="magma", \
               extent=extents,origin="lower",vmin=v2n,vmax=v2p)
    if (ilev==0):
      cbaxes = fig.add_axes([0.40, 0.8, 0.23, 0.02])
      cb=fig.colorbar(im2, ax=ax2,orientation="horizontal", cax=cbaxes)
      cb.ax.tick_params(labelsize=8)
      tick_locator = ticker.MaxNLocator(nbins=4)
      cb.locator = tick_locator
      cb.update_ticks()

    if ilev==0:
      bn=bp=Pmag.min(); bp=Pmag.max()
    im3 = ax3.imshow(Pmag,interpolation="nearest",cmap="plasma", \
               extent=extents,origin="lower",vmin=bn,vmax=bp)
    if (ilev==0):
      cbaxes = fig.add_axes([0.67, 0.8, 0.23, 0.02]) 
      cb=fig.colorbar(im3, ax=ax3,orientation="horizontal", cax=cbaxes,pad=0.0)
      cb.ax.tick_params(labelsize=8)
      tick_locator = ticker.MaxNLocator(nbins=4)
      cb.locator = tick_locator
      cb.update_ticks()

  dmax=str("%.2f" % (d.max()))
  dmin=str("%.2f" % (d.min()))
  tm = str("%.1f" % t)
  st = "(a) Density, $\\rho \in$[" +dmin + "," + dmax +"]"
  ax1.text(-0.6,0.775,st,color="black",fontsize=8)

  V2 = 0.5*np.sqrt(vx[0]*vx[0] + vy[0]*vy[0] + vz[0]*vz[0])
  vmax=str("%.2f" % (V2.max()))
  vmin=str("%.2f" % (V2.min()))
  st = "(b) Velocity, $\left| \mathbf{v} \\right| \in$[" +vmin + "," + vmax +"]"
  ax2.text(-0.6,0.775,st,color="black",fontsize=8)

  Pm = (bx[0]*bx[0] + by[0]*by[0] + bz[0]*bz[0])/(8.0*np.pi)
  bmax=str("%.2f" % (Pm.max()))
  bmin=str("%.2f" % (Pm.min()))
  st = "(c) $P_\mathrm{mag}\equiv\left| B \\right|^2/2 \in$[" +bmin + "," + bmax +"]"
  ax3.text(-0.6,0.775,st,color="black",fontsize=8)
  
  fn = file[0][0:12] + "_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".png"
  plt.savefig(fn,bbox_inches="tight",dpi=300)

quit()
