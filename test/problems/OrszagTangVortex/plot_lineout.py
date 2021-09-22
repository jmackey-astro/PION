#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/jm/code/pypion/silo/lib")
sys.path.append("/home/jm/code/pypion/Library")
import Silo
import Plotting_Classes as ppion
import glob

plt.rcParams["font.weight"] = "normal"
plt.rcParams["font.size"] = "12"
plt.rcParams["font.family"] = "serif"
plt.rcParams['lines.linewidth'] = 1
plt.rcParams["mathtext.fontset"] = "cm"

# plot 3 different resolutions at t=0.5
f0=sorted(glob.glob("OrszagTang_n128_b3.33m1.0_devel_0000.*.silo"))
f1=sorted(glob.glob("OrszagTang_n256_b3.33m1.0_devel_0000.*.silo"))
f2=sorted(glob.glob("OrszagTang_n512_b3.33m1.0_devel_0000.*.silo"))
f = [f0[5],f1[5],f2[5]]
# two plots, lineouts at y=0.3125 and 0.4277
fig, (ax1,ax2) = plt.subplots(2,1)
styles=["r:","b--","k-"]
ct=0

# loop over files to plot each one in turn.
for j in f:
  print(j)
  file = [j]
  print(file)
  
  dataio=ppion.Plotting2d(file)
  n = dataio.nlevels()
  c = dataio.cycle()
  level_ng  = dataio.ngrid()
  P = dataio.get_2Darray("Pressure")
  p = P['data']
  t = P['sim_time']
  print("time=",t)
  xmax = P['max_extents']
  xmin = P['min_extents']

  ilev=0
  #print("ilev=",ilev)
  level_min = xmin[ilev]
  level_max = xmax[ilev]
  dy=((level_max[1]-level_min[1]))/level_ng[1]
  dx=((level_max[0]-level_min[0]))/level_ng[0]
  x = np.arange(level_min[0]+0.5*dx,level_max[0]+0.5*dx,dx)
  #print(x)
  p = np.array(p)
  pnorm = 4*np.pi

  # Get a lineout along x for y=0.3125 (as near as possible), save as p1
  yseek=0.3125
  iseek=0
  diff=1e99
  while (level_min[1]+(iseek+0.5)*dy <yseek):
    iseek = iseek+1
  #print(level_min[1]+(iseek-0.5)*dy- yseek,level_min[1]+(iseek+0.5)*dy- yseek)
  #print("iseek=",iseek,", y=",level_min[1]+(iseek+0.5)*dy)

  # take a weighted average of the two rows nearest the y-value we seek:
  offup=level_min[1]+(iseek+0.5)*dy - yseek
  offdn=yseek - level_min[1]+(iseek+0.5)*dy
  p1 = (offup*p[ilev][iseek,:] + offdn*p[ilev][iseek-1,:])/(offup+offdn)/pnorm
  p1=p[ilev][iseek,:]/pnorm
  #print("p1",p1)

  # Get a lineout along x for y=0.4277 (as near as possible), save as p1
  yseek=0.4277
  iseek=0
  diff=1e99
  while (level_min[1]+(iseek+0.5)*dy <yseek):
    iseek = iseek+1
  #print(level_min[1]+(iseek-0.5)*dy- yseek,level_min[1]+(iseek+0.5)*dy- yseek)
  #print("iseek=",iseek,", y=",level_min[1]+(iseek+0.5)*dy)
  
  # take a weighted average of the two rows nearest the y-value we seek:
  offup=level_min[1]+(iseek+0.5)*dy - yseek
  offdn=yseek - level_min[1]+(iseek+0.5)*dy
  p2 = (offup*p[ilev][iseek,:] + offdn*p[ilev][iseek-1,:])/(offup+offdn)/pnorm
  p2=p[ilev][iseek,:]/pnorm
  #print("p2",p2)

  ax1.set_xlim(level_min[0], level_max[0])
  ax1.set_ylim(0,0.3)
  ax1.set_xlabel("$x$")
  ax1.set_ylabel("$p$")
  ax1.grid()
  #ax1.xaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
  #ax1.tick_params(labelsize=8)
  
  ax2.set_xlim(level_min[0], level_max[0])
  ax2.set_ylim(0,0.5)
  ax2.set_xlabel("$x$")
  ax2.set_ylabel("$p$")
  ax2.grid()
  #ax2.yaxis.set_ticklabels([])
  #ax2.xaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
  #ax2.tick_params(labelsize=8)

  ax1.plot(x,p1,styles[ct])
  ax2.plot(x,p2,styles[ct])
  ct = ct+1


tm = str("%.2f" % t)
st = "(a) Pressure $y=0.3125$, $t=$"+tm
ax1.text(0.02,0.26,st,color="black",fontsize=10)
st = "(b) Pressure $y=0.4277$, $t=$"+tm
ax2.text(0.02,0.45,st,color="black",fontsize=10)
fn = file[0][0:10] + "_l"+str(n)+"n"+str(level_ng[0])+"c"+str(c)+".png"
plt.savefig(fn,bbox_inches="tight",dpi=300)


quit()
