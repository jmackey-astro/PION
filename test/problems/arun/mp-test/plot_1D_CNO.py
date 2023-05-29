# -*- coding: iso-8859-15 -*-

import sys
sys.path.append("/home/tony/Desktop/GitLab_PION/pypion/silo/lib")
sys.path.append("/home/tony/Desktop/GitLab_PION/pypion/Library")
import Plotting_Classes as pypion
from ReadData import ReadData
from SiloHeader_data import OpenData

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

for i in range(0,len(files)):
  print(i,files[i])
  #dataio=OpenData([files[i]])
  dataio = ReadData([files[i]])
  n = dataio.nlevels()
  c = dataio.cycle()
  #print(n,c)
  D = dataio.get_1Darray("Density")
  ro = D['data'][0]
  time = (D['sim_time'] * u.s).to(u.Myr)
  print("time=",time)
  xmax = (D['max_extents'] * u.cm).to(u.pc)
  xmin = (D['min_extents'] * u.cm).to(u.pc)
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

  pg = dataio.get_1Darray("Pressure")['data'][0]
  vx = dataio.get_1Darray("VelocityX")['data'][0]
  vy = dataio.get_1Darray("VelocityY")['data'][0]
  vz = dataio.get_1Darray("VelocityZ")['data'][0]
  tm = dataio.get_1Darray("Temperature")['data'][0]

  # Read element abundances
  xh  = dataio.get_1Darray("Tr000_X_H")['data'][0]
  xhe = dataio.get_1Darray("Tr001_X_He")['data'][0]
  xc  = dataio.get_1Darray("Tr002_X_C")['data'][0]
  xn  = dataio.get_1Darray("Tr003_X_N")['data'][0]
  xo  = dataio.get_1Darray("Tr004_X_O")['data'][0]

  # Read ion abundances
  hp   = dataio.get_1Darray("Tr005_H1p")['data'][0]
  he1p = dataio.get_1Darray("Tr006_He1p")['data'][0]
  he2p = dataio.get_1Darray("Tr007_He2p")['data'][0]
  c1p  = dataio.get_1Darray("Tr008_C1p")['data'][0]
  c2p  = dataio.get_1Darray("Tr009_C2p")['data'][0]
  c3p  = dataio.get_1Darray("Tr010_C3p")['data'][0]
  c4p  = dataio.get_1Darray("Tr011_C4p")['data'][0]
  c5p  = dataio.get_1Darray("Tr012_C5p")['data'][0]
  c6p  = dataio.get_1Darray("Tr013_C6p")['data'][0]
  n1p  = dataio.get_1Darray("Tr014_N1p")['data'][0]
  n2p  = dataio.get_1Darray("Tr015_N2p")['data'][0]
  n3p  = dataio.get_1Darray("Tr016_N3p")['data'][0]
  n4p  = dataio.get_1Darray("Tr017_N4p")['data'][0]
  n5p  = dataio.get_1Darray("Tr018_N5p")['data'][0]
  n6p  = dataio.get_1Darray("Tr019_N6p")['data'][0]
  n7p  = dataio.get_1Darray("Tr020_N7p")['data'][0]  
  o1p  = dataio.get_1Darray("Tr021_O1p")['data'][0]
  o2p  = dataio.get_1Darray("Tr022_O2p")['data'][0]
  o3p  = dataio.get_1Darray("Tr023_O3p")['data'][0]
  o4p  = dataio.get_1Darray("Tr024_O4p")['data'][0]
  o5p  = dataio.get_1Darray("Tr025_O5p")['data'][0]
  o6p  = dataio.get_1Darray("Tr026_O6p")['data'][0]
  o7p  = dataio.get_1Darray("Tr027_O7p")['data'][0]
  o8p  = dataio.get_1Darray("Tr028_O8p")['data'][0]

  # Calculate neutral species abundances from conservation equation
  h0  = xh  - hp
  he0 = xhe -(he1p+he2p)
  c0  = xc  -(c1p+c2p+c3p+c4p+c5p+c6p)
  n0  = xn  -(n1p+n2p+n3p+n4p+n5p+n6p+n7p)
  o0  = xn  -(o1p+o2p+o3p+o4p+o5p+o6p+o7p+o8p)
  np.clip(h0,  1.0e-10, 1.0, out=h0)
  np.clip(he0, 1.0e-10, 1.0, out=he0)
  np.clip(c0,  1.0e-10, 1.0, out=c0)
  np.clip(n0,  1.0e-10, 1.0, out=n0)
  np.clip(o0,  1.0e-10, 1.0, out=o0)

# convert radius to pc, velocity to km/s
#  x  = x*((1.0*u.cm).to(u.pc)).value
  vr = vx*1.0e-5
  nh = ro*xh/(apc.m_p.cgs.value)

# Plot H / He ions
  fig = plt.figure()
  plt.plot(np.log10(tm),np.log10(h0),"r-",label="$y(\mathrm{H}^{0})$")
  plt.plot(np.log10(tm),np.log10(hp),"r--",label="$y(\mathrm{H}^{+})$")
  plt.plot(np.log10(tm),np.log10(he0),"b-",label="$y(\mathrm{He}^{0})$")
  plt.plot(np.log10(tm),np.log10(he1p),"b--",label="$y(\mathrm{He}^{+})$")
  plt.plot(np.log10(tm),np.log10(he2p),"b-.",label="$y(\mathrm{He}^{2+})$")
  plt.ylim(-13,0)

  plt.tick_params(labelsize=16)
  plt.grid()
  #plt.legend(fontsize=12, loc=([0.25,0.25]))
  plt.legend(fontsize=12, loc="lower right")
  plt.ylabel("$\log_{10}$ quantities",fontsize=16)
  plt.xlabel("$\log_{10}T$ (K)",fontsize=16)
  s = "$t=$" + f"{time:0.04f}"
  plt.text(1.5,1,s,color="black",fontsize=14)
#plt.savefig("Gal_M30V400_d1l1n1024.pdf", bbox_inches="tight", density=600)
  iy = str(i).zfill(5)
  opf = img_path+"/H&He_"+fbase+"."+iy+".png"
  plt.savefig(opf, bbox_inches="tight")
  plt.close(fig)
  del fig

# Plot Carbon ions
  fig = plt.figure()
  plt.plot(np.log10(tm),np.log10(c0),"b-",label="$y(\mathrm{C}^{0})$")
  plt.plot(np.log10(tm),np.log10(c1p),"b--",label="$y(\mathrm{C}^{+})$")
  plt.plot(np.log10(tm),np.log10(c2p),"r-.",label="$y(\mathrm{C}^{2+})$")
  plt.plot(np.log10(tm),np.log10(c3p),"b:",label="$y(\mathrm{C}^{3+})$")
  plt.plot(np.log10(tm),np.log10(c4p),"g-",label="$y(\mathrm{C}^{4+})$")
  plt.plot(np.log10(tm),np.log10(c5p),"g--",label="$y(\mathrm{C}^{5+})$")
  plt.plot(np.log10(tm),np.log10(c6p),"g-.",label="$y(\mathrm{C}^{6+})$")
  #plt.xlim(0,70)
  plt.ylim(-13,0)

  plt.tick_params(labelsize=16)
  plt.grid()
  #plt.legend(fontsize=12, loc=([0.25,0.25]))
  plt.legend(fontsize=12, loc="lower right")
  plt.ylabel("$\log_{10}$ quantities",fontsize=16)
  plt.xlabel("$\log_{10}T$ (K)",fontsize=16)
  s = "$t=$" + f"{time:0.04f}"
  plt.text(1.5,1,s,color="black",fontsize=14)
#plt.savefig("Gal_M30V400_d1l1n1024.pdf", bbox_inches="tight", density=600)
  iy = str(i).zfill(5)
  opf = img_path+"/C_"+fbase+"."+iy+".png"
  plt.savefig(opf, bbox_inches="tight")
  plt.close(fig)
  del fig

# Plot Nitrogen Iona
  fig = plt.figure()
  plt.plot(np.log10(tm),np.log10(n0),"b-",label="$y(\mathrm{N}^{0})$")
  plt.plot(np.log10(tm),np.log10(n1p),"r-",label="$y(\mathrm{N}^{+})$")
  plt.plot(np.log10(tm),np.log10(n2p),"b-.",label="$y(\mathrm{N}^{2+})$")
  plt.plot(np.log10(tm),np.log10(n3p),"b:",label="$y(\mathrm{N}^{3+})$")
  plt.plot(np.log10(tm),np.log10(n4p),"g-",label="$y(\mathrm{N}^{4+})$")
  plt.plot(np.log10(tm),np.log10(n5p),"g--",label="$y(\mathrm{N}^{5+})$")
  plt.plot(np.log10(tm),np.log10(n6p),"g-.",label="$y(\mathrm{N}^{6+})$")
  plt.plot(np.log10(tm),np.log10(n7p),"g:",label="$y(\mathrm{N}^{7+})$")
  plt.ylim(-13,0)

  plt.tick_params(labelsize=16)
  plt.grid()
  plt.legend(fontsize=12, loc="lower right")
  plt.ylabel("$\log_{10}$ quantities",fontsize=16)
  plt.xlabel("$\log_{10}T$ (K)",fontsize=16)
  s = "$t=$" + f"{time:0.04f}"
  plt.text(1.5,1,s,color="black",fontsize=14)
  iy = str(i).zfill(5)
  opf = img_path+"/N_"+fbase+"."+iy+".png"
  plt.savefig(opf, bbox_inches="tight")
  plt.close(fig)
  del fig
  
  # Plot Oxygen Iona
  fig = plt.figure()
  plt.plot(np.log10(tm),np.log10(o0),"b-",label="$y(\mathrm{O}^{0})$")
  plt.plot(np.log10(tm),np.log10(o1p),"r-",label="$y(\mathrm{O}^{+})$")
  plt.plot(np.log10(tm),np.log10(o2p),"g-",label="$y(\mathrm{O}^{2+})$")
  plt.plot(np.log10(tm),np.log10(o3p),"b--",label="$y(\mathrm{O}^{3+})$")
  plt.plot(np.log10(tm),np.log10(o4p),"b-.",label="$y(\mathrm{O}^{4+})$")
  plt.plot(np.log10(tm),np.log10(o5p),"b:",label="$y(\mathrm{O}^{5+})$")
  plt.plot(np.log10(tm),np.log10(o6p),"g-",label="$y(\mathrm{O}^{6+})$")
  plt.plot(np.log10(tm),np.log10(o7p),"g--",label="$y(\mathrm{O}^{7+})$")
  plt.plot(np.log10(tm),np.log10(o8p),"g-.",label="$y(\mathrm{O}^{8+})$")
  plt.ylim(-13,0)

  plt.tick_params(labelsize=16)
  plt.grid()
  plt.legend(fontsize=12, loc="lower right")
  plt.ylabel("$\log_{10}$ quantities",fontsize=16)
  plt.xlabel("$\log_{10}T$ (K)",fontsize=16)
  s = "$t=$" + f"{time:0.04f}"
  plt.text(1.5,1,s,color="black",fontsize=14)
  iy = str(i).zfill(5)
  opf = img_path+"/O_"+fbase+"."+iy+".png"
  plt.savefig(opf, bbox_inches="tight")
  plt.close(fig)
  del fig


  dataio.close()
  del dataio
  del x, ro, nh, pg, vx, vr, vy, vz, tm, xh, xhe, hp
  del xmin, xmax, dx, ng, time, D, c, n

quit()



