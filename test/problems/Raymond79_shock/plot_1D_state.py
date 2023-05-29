# -*- coding: iso-8859-15 -*-

import sys
sys.path.insert(0,"/home/jmackey/.local/silo/lib")
#sys.path.insert(0,"/home/jm/code/pypion/Library")
#sys.path.insert(0,"/home/jmackey/.local/silo/lib")
#sys.path.insert(0,"/mnt/local/jm/pion_python/src/pypion/")
#sys.path.insert(0,"/home/jmackey/code/pypion/silo/lib")
#sys.path.insert(0,"/home/jmackey/code/pypion/Library")
import Silo
#import Plotting_Classes as ppion
#sys.path.append("/mnt/massive-stars/share/pypion/silo/lib")
#sys.path.append("/mnt/massive-stars/share/pypion/Library")
#import Plotting_Classes as pypion
from pypion.ReadData import ReadData
from pypion.SiloHeader_data import OpenData

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('text', usetex=True)
#plt.rc('font',**{'size': 20})
#plt.rc('lines', linewidth=2)
import argparse
from os import listdir
import glob

import astropy.units as u
from astropy import constants as apc


plt.rcParams["font.weight"] = "normal"
from matplotlib.font_manager import FontProperties

font = FontProperties()
font.set_family('sans-serif')
font.set_name('stixsans')
font.set_style('italic')
font.set_weight('light')
font.set_size(12)

parser = argparse.ArgumentParser()
parser.add_argument("path", type=str)
parser.add_argument("fbase", type=str)
parser.add_argument("img_path", type=str)
args = parser.parse_args()
path = args.path
fbase = args.fbase
img_path = args.img_path

files = sorted(glob.glob(path+"/"+fbase+".*.silo"))

for i in range(len(files)):
  print(i,files[i])
  #dataio=OpenData([files[i]])
  dataio = ReadData([files[i]])
  n = dataio.nlevels()
  c = dataio.cycle()
  #print(n,c)
  D = dataio.get_1Darray("Density")
  ro = D['data'][0]
  time = (D['sim_time'] * u.s).to(u.kyr)
  print("time=",time)
  xmax = (D['max_extents'] * u.cm).to(u.AU)
  xmin = (D['min_extents'] * u.cm).to(u.AU)
  ng  = dataio.ngrid()
  dx = (xmax-xmin)/ng
  xmin = xmin[0]
  xmax = xmax[0]
  ng = ng[0]
  dx = dx[0]
  #print(xmin[0],xmax[0],dx[0])
  x0 = xmin[0]+0.5*dx[0]
  xn = xmax[0]-0.5*dx[0]
  x = np.linspace(x0,xn,ng)
  #print(xmin,xmax,ng,dx,x)
  #print(x)

  #pg = dataio.get_1Darray("Pressure")['data'][0]
  vx = dataio.get_1Darray("VelocityX")['data'][0]
  #vy = dataio.get_1Darray("VelocityY")['data'][0]
  #vz = dataio.get_1Darray("VelocityZ")['data'][0]
  tm = dataio.get_1Darray("Temperature")['data'][0]
  xh = dataio.get_1Darray("Tr000_X_H")['data'][0]
  #xhe = dataio.get_1Darray("Tr001_X_He")['data'][0]
  hp = dataio.get_1Darray("Tr005_H1p")['data'][0]
  hep = dataio.get_1Darray("Tr006_He1p")['data'][0]
  he2p = dataio.get_1Darray("Tr007_He2p")['data'][0]
  #ct = dataio.get_1Darray("Tr008_WIND")['data'][0]
  
# convert radius to pc, velocity to km/s
#  x  = x*((1.0*u.cm).to(u.pc)).value
  vr = vx*1.0e-5
  nh = ro*xh/(apc.m_p.cgs.value)


  fig = plt.figure()
#plt.loglog(x,nh,"k-",label="$n_\mathrm{H}\, \left(\mathrm{cm}^{-3}\\right)$")
#plt.loglog(x,np.fabs(vr),"b--",label="$\left|v_r\\right|\, \left(\mathrm{km\,s}^{-1}\\right)$")
#plt.loglog(x,tm*1.0e-5,"r-.",label="$T\, (10^5\,\mathrm{K})$")
#plt.loglog(x,hp,"k:",label="$y(\mathrm{H}^{+})$")
#plt.xlim(0.1,50)
#plt.ylim(1e-5,0.9e4)

  plt.plot(x,np.log10(nh),"k-",label="$n_\mathrm{H}\, \left(\mathrm{cm}^{-3}\\right)$")
  plt.plot(x,np.log10(np.fabs(vr)),"b--",label="$\left|v_r\\right|\, \left(\mathrm{km\,s}^{-1}\\right)$")
  plt.plot(x,np.log10(tm*1.0e-5),"r-.",label="$T\, (10^5\,\mathrm{K})$")
  plt.plot(x,np.log10(hp),"k:",label="$y(\mathrm{H}^{+})$")
  plt.plot(x,np.log10(hep),"m--",label="$y(\mathrm{He}^{+})$")
  plt.plot(x,np.log10(he2p),"m:",label="$y(\mathrm{He}^{2+})$")
  #plt.xlim(0,5)
  plt.ylim(-5.5,4)

  plt.tick_params(labelsize=16)
  plt.grid()
  #plt.legend(fontsize=12, loc=([0.25,0.25]))
  plt.legend(fontsize=12, loc="lower right")
  plt.ylabel("various units",fontsize=16)
  plt.xlabel("$x$ (AU)",fontsize=16)
  s = "$t=$" + f"{time:0.03f}"
  #plt.text(42,3.4,s,color="black",fontsize=14)
  #plt.xlim(0,5)
  plt.text(42,3.4,s,color="black",fontsize=14)
#plt.savefig("Gal_M30V400_d1l1n1024.pdf", bbox_inches="tight", density=600)
  iy = str(i).zfill(5)
  opf = img_path+"/"+fbase+"."+iy+".png"
  plt.savefig(opf, bbox_inches="tight")
  plt.close(fig)
  del fig
  dataio.close()
  del dataio
  #del x, ro, nh, pg, vx, vr, vy, vz, tm, xh, xhe, hp
  #del xmin, xmax, dx, ng, time, D, c, n

quit()



