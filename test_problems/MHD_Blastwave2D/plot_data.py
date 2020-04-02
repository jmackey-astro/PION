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


#plt.rc('font',**{'family':'sans-serif','sans-serif':['Bitstream Vera Sans']})
#plt.rc('text', usetex=True)
#plt.rc('font',**{'size': 12})
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
font.set_size(8)

solver="HLLD_NG"
#solver="HLLD_UG"
if solver=="HLL":
  files=[""]
elif solver=="HLLD_NG":
  #files=glob.glob("HLLD_B001_n040_l2_level00_0000.*.silo")
  files=glob.glob("HLLD_B010_n040_l2_level00_0000.*.silo")
  #files=glob.glob("HLLD_B001_n200_l2_level00_0000.*.silo")
  files=sorted(glob.glob("BW2d_StoneMHD_B010_n200_level00.000*.silo"))
elif solver=="HLLD_UG":
  files=glob.glob("HLLD_n200_B1_ug_0000.*.silo")
elif solver=="RCV":
  files==glob.glob("RCV*_level00_0000.*.silo")
else:
  files=[""]
  #files=["FieldLoop032_0000.00000000.silo", "FieldLoop032_0000.00000132.silo", "FieldLoop032_0000.00000264.silo"]

for file in files:
  print file
  fig, (ax1,ax2,ax3) = plt.subplots(1,3)
  
  dataio=ppion.Plotting2d(file)
  n = dataio.nlevels()
  c = dataio.cycle()
  t = dataio.sim_time().value
  level_ng  = dataio.ngrid()

  dn=0.0; dp=0.0;
  v2n=0.0; v2p=0.0;
  bn=0.0; bp=0.0;
  
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
    d = dataio.reshaped_parameter2d("Density")
    #d = dataio.reshaped_parameter2d("Pressure")
    vx = dataio.reshaped_parameter2d("VelocityX")
    vy = dataio.reshaped_parameter2d("VelocityY")
    vz = dataio.reshaped_parameter2d("VelocityZ")
    bx = dataio.reshaped_parameter2d("MagneticFieldX")
    by = dataio.reshaped_parameter2d("MagneticFieldY")
    bz = dataio.reshaped_parameter2d("MagneticFieldZ")
    B = np.sqrt(4.0*np.pi*(bx*bx + by*by + bz*bz))
    V2 = 0.5*np.sqrt(vx*vx + vy*vy + vz*vz)
    Pmag=B**2/(8.0*np.pi)
    
    level_min = dataio.level_min()
    level_max = dataio.level_max()
    dy=((level_max[1]-level_min[1])).value/level_ng[1]
    dx=((level_max[0]-level_min[0])).value/level_ng[0]
    y, x = np.mgrid[slice(level_min[1].value+0.5*dy, level_max[1].value+0.5*dy, dy),
                    slice(level_min[0].value+0.5*dx, level_max[0].value+0.5*dx, dx)]
    if (ilev==0):
      ax1.set_xlim(level_min[0].value, level_max[0].value)
      ax1.set_ylim(level_min[1].value, level_max[1].value)
      #ax1.set_xlim(-1.0,1.0)
      #ax1.set_ylim(-0.5,0.5)
      ax1.set_xlabel("x",fontsize=10, fontweight='bold')
      ax1.set_ylabel("y",fontsize=10, fontweight='bold')
      ax1.xaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
      ax1.tick_params(labelsize=8)
      
      ax2.set_xlim(level_min[0].value, level_max[0].value)
      ax2.set_ylim(level_min[1].value, level_max[1].value)
      #ax2.set_xlim(-1.0,1.0)
      #ax2.set_ylim(-0.5,0.5)
      ax2.set_xlabel("x",fontsize=10, fontweight='bold')
      #ax2.set_ylabel("y",fontsize=12, fontweight='bold')
      ax2.yaxis.set_ticklabels([])
      ax2.xaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
      ax2.tick_params(labelsize=8)
      
      ax3.set_xlim(level_min[0].value, level_max[0].value)
      ax3.set_ylim(level_min[1].value, level_max[1].value)
      #ax3.set_xlim(-1.0,1.0)
      #ax3.set_ylim(-0.5,0.5)
      ax3.set_xlabel("x",fontsize=10, fontweight='bold')
      #ax3.set_ylabel("y",fontsize=12, fontweight='bold')
      ax3.yaxis.set_ticklabels([])
      ax3.xaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
      ax3.tick_params(labelsize=8)
      
    #if (ilev>0) or (n==1):
      #plt.pcolormesh(x, y, d, cmap="viridis",linewidth=0,rasterized=True)
      clev = MaxNLocator(nbins=21).tick_values(d.min(), d.max())
      #plt.contour(x, y, B, colors="white",levels=clev,linewidths=0.5,aspect=1)
      #clev = [0.0,0.002,0.004,0.006,0.008,0.010,0.012,0.014]
      print clev
      ct1 = ax1.contour(x, y, d, colors="white",levels=clev,linewidths=0.25)
      ct1.monochrome = True

      clev = MaxNLocator(nbins=7).tick_values(V2.min(), V2.max())
      #clev = [-0.016 -0.008  0.     0.008  0.016  0.024  0.032  0.04 ]
      print clev
      ct2 = ax2.contour(x, y, V2, colors="white",levels=clev,linewidths=0.5)
      ct2.monochrome = True

      clev = MaxNLocator(nbins=7).tick_values(B.min(), B.max())
      #clev = [-0.016 -0.008  0.     0.008  0.016  0.024  0.032  0.04 ]
      print clev
      ct3 = ax3.contour(x, y, B, colors="white",levels=clev,linewidths=0.5)
      ct3.monochrome = True

    extents = [level_min[0].value, level_max[0].value, \
               level_min[1].value,level_max[1].value]

    #clev = MaxNLocator(nbins=11).tick_values(d.min(),d.max())
    #print(clev)
    #im1 =ax1.contourf(x, y, d, cmap="viridis",levels=clev)
    if ilev==0:
      dn=0.0; dp=d.max()
    im1 = ax1.imshow(d,interpolation="nearest",cmap="viridis", \
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
    print "lims v2 plot lev=",ilev,v2n,v2p,V2.min(),V2.max()
    im2 = ax2.imshow(V2,interpolation="nearest",cmap="magma", \
               extent=extents,origin="lower",vmin=v2n,vmax=v2p)
    if (ilev==0):
      #fig.colorbar(im2, ax=ax2, shrink=0.6)
      cbaxes = fig.add_axes([0.40, 0.8, 0.23, 0.02])
      #cbar_kws=dict(ticks=(0.0,0.08,0.16,0.24,0.32))
      #cb=fig.colorbar(im2, ax=ax2,orientation="horizontal", cax=cbaxes,**dict(cbar_kws))
      cb=fig.colorbar(im2, ax=ax2,orientation="horizontal", cax=cbaxes)
      cb.ax.tick_params(labelsize=8)
      #cb.ax.ticks([0.0,0.08,0.16,0.24,0.32])
      tick_locator = ticker.MaxNLocator(nbins=4)
      cb.locator = tick_locator
      cb.update_ticks()

    if ilev==0:
      bn=0.0; bp=B.max()
    im3 = ax3.imshow(B,interpolation="nearest",cmap="plasma", \
               extent=extents,origin="lower",vmin=bn,vmax=bp)
    if (ilev==0):
      #fig.colorbar(im3, ax=ax3, shrink=0.4)
      cbaxes = fig.add_axes([0.67, 0.8, 0.23, 0.02]) 
      cb=fig.colorbar(im3, ax=ax3,orientation="horizontal", cax=cbaxes,pad=0.0)
      cb.ax.tick_params(labelsize=8)
      tick_locator = ticker.MaxNLocator(nbins=4)
      cb.locator = tick_locator
      cb.update_ticks()

    dataio.close()
    del dataio


  dmax=str("%.2f" % (d.max()))
  dmin=str("%.2f" % (d.min()))
  tm = str("%.3f" % t)
  #st = "$P_\mathrm{mag}$, max=" +dmax + " min=" + dmin
  st = "(a) $\\rho$, t=" + tm + "\nmax=" +dmax + " min=" + dmin
  #st = "Density, t="+tm
  #ax1.text(-0.75,0.42,st,color="white",fontsize=10)
  #if c==0:
  ax1.text(-0.475,0.55,st,color="white",fontsize=8,weight='bold')
  #else:
  #  ax1.text(-0.95,0.42,"(c)",color="white",fontsize=10,weight='bold')

  st = "$\left| B \\right|$, max=" +dmax + " min=" + dmin

  dmax=str("%.4f" % (V2.max()))
  dmin=str("%.4f" % (V2.min()))
  tm = str("%.5f" % t)
  st = "$\\nabla\\times B$, max=" +dmax + " min=" + dmin
  #st = "Density, t="+tm
  #ax2.text(-0.75,0.42,st,color="black",fontsize=10)
  #if c==0:
  #  ax2.text(-0.95,0.42,"(b)",color="black",fontsize=10,weight='bold')
  #else:
  #  ax2.text(-0.95,0.42,"(d)",color="black",fontsize=10,weight='bold')

  #for ax in fig.get_axes():
  #  ax.label_outer()

  #plt.show()
  fn = file[0:26] + "_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".png"
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
