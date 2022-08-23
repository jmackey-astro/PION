# -*- coding: iso-8859-15 -*-

import sys
#sys.path.insert(0,"/home/jm/code/pypion/silo/lib")
sys.path.insert(0,"/home/jm/.local/silo/lib")
import Silo
#sys.path.insert(0,"/home/jm/code/pypion/Library")
from pypion.ReadData import ReadData

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
# ignore divide-by-zero warnings :)
np.seterr(divide='ignore')

parser = argparse.ArgumentParser()
parser.add_argument("path", type=str)
parser.add_argument("fbase", type=str)
parser.add_argument("img_path", type=str)
args = parser.parse_args()
path = args.path
sim = args.fbase
img_path = args.img_path

# Poor way of checking how many grid levels there are: I always
# name files with number of dimensions "d" and levels "l"
if sim.find("d1l2") > 0:
  files0 = sorted(glob.glob(path+"/"+sim+"_level00_0000.*.silo"))
  files1 = sorted(glob.glob(path+"/"+sim+"_level01_0000.*.silo"))
  lev=2
elif sim.find("d1l1") > 0:
  files0 = sorted(glob.glob(path+"/"+sim+"_0000.*.silo"))
  lev=1
elif sim.find("d1l3") > 0:
  files0 = sorted(glob.glob(path+"/"+sim+"_level00_0000.*.silo"))
  files1 = sorted(glob.glob(path+"/"+sim+"_level01_0000.*.silo"))
  files2 = sorted(glob.glob(path+"/"+sim+"_level02_0000.*.silo"))
  lev=3
elif sim.find("d1l4") > 0:
  files0 = sorted(glob.glob(path+"/"+sim+"_level00_0000.*.silo"))
  files1 = sorted(glob.glob(path+"/"+sim+"_level01_0000.*.silo"))
  files2 = sorted(glob.glob(path+"/"+sim+"_level02_0000.*.silo"))
  files3 = sorted(glob.glob(path+"/"+sim+"_level03_0000.*.silo"))
  lev=4
elif sim.find("d1l5") > 0:
  files0 = sorted(glob.glob(path+"/"+sim+"_level00_0000.*.silo"))
  files1 = sorted(glob.glob(path+"/"+sim+"_level01_0000.*.silo"))
  files2 = sorted(glob.glob(path+"/"+sim+"_level02_0000.*.silo"))
  files3 = sorted(glob.glob(path+"/"+sim+"_level03_0000.*.silo"))
  files4 = sorted(glob.glob(path+"/"+sim+"_level04_0000.*.silo"))
  lev=5
else:
  files0 = sorted(glob.glob(path+"/"+sim+"_0000.*.silo"))
  lev=1

if (len(files0)<1):
  print("can't find requested files")
  quit()


for i in range(len(files0)):
  print(i,files0[i])
  # PyPion needs a list of all files to read:
  if lev==1:
    files = [files0[i]]
  elif lev==2:
    files =[files0[i],files1[i]]
  elif lev==3:
    files = [files0[i],files1[i],files2[i]]
  elif lev==4:
    files = [files0[i],files1[i],files2[i],files3[i]]
  elif lev==5:
    files = [files0[i],files1[i],files2[i],files3[i],files4[i]]
  else:
    print("bad number of levels")
    quit()
  dataio = ReadData(files)

  n = dataio.nlevels()
  c = dataio.cycle()
  #print(n,c)
  D = dataio.get_1Darray("Density")
  vP = dataio.get_1Darray("Pressure")['data']
  vV = dataio.get_1Darray("VelocityX")['data']
  vE = dataio.get_1Darray("Temperature")['data']
#  xh = dataio.get_1Darray("Tr000_X_H")['data'][0]
#  xhe = dataio.get_1Darray("Tr001_X_He")['data'][0]
  vhp = dataio.get_1Darray("Tr000_H1p")['data']
  vT = dataio.get_1Darray("Tr001_WIND")['data']
  # NG_Mask = 0 if cell is not a leaf, =1 if cell is a leaf.
  if lev>1:
    vM = dataio.get_1Darray("NG_Mask")['data']
  vD = D['data']  # this is an array of arrays, one for each grid level

  time = (D['sim_time'] * u.s).to(u.kyr)
  print("time=",time)
  xmax = (D['max_extents'] * u.cm).to(u.pc)
  xmin = (D['min_extents'] * u.cm).to(u.pc)
  ng  = dataio.ngrid()

  fig = plt.figure()
  plt.ylabel("various units",fontsize=16)
  plt.xlabel("$x$ (pc)",fontsize=16)
  plt.xlim(0, xmax[0][0].value)
  plt.ylim(-4,8.5)
  s = "$t=$" + f"{time:0.03f}"
  plt.text(0.02,4.4,s,color="black",fontsize=14)
  plt.tick_params(labelsize=16)
  plt.grid()

  # Loop over levels
  for ilev in range(n):
    lmin = xmin[ilev].value
    lmax = xmax[ilev].value
    dx = (lmax[0]-lmin[0])/ng[0]
    x0 = lmin[0]+0.5*dx
    xn = lmax[0]-0.5*dx
    x = np.linspace(x0,xn,ng[0])
    
    if lev>1:
      ro = vD[ilev]*vM[ilev]  # multiply by mask to zero non-leaf cells
      pg = vP[ilev]*vM[ilev]
      vx = vV[ilev]*vM[ilev]
      tt = vE[ilev]*vM[ilev]
      tr = vT[ilev]*vM[ilev]
      hp = vhp[ilev]*vM[ilev]
    else:
      ro = vD[ilev]
      pg = vP[ilev]
      vx = vV[ilev]
      tt = vE[ilev]
      tr = vT[ilev]
      hp = vhp[ilev]

    xh = 0.72
    xhe = 0.26
    vr = vx*1.0e-5  # convert to km/s
    nh = ro*xh/(apc.m_p.cgs.value)  # convert density to n(H)


#plt.loglog(x,nh,"k-",label="$n_\mathrm{H}\, \left(\mathrm{cm}^{-3}\\right)$")
#plt.loglog(x,np.fabs(vr),"b--",label="$\left|v_r\\right|\, \left(\mathrm{km\,s}^{-1}\\right)$")
#plt.loglog(x,tm*1.0e-5,"r-.",label="$T\, (10^5\,\mathrm{K})$")
#plt.loglog(x,hp,"k:",label="$y(\mathrm{H}^{+})$")
#plt.xlim(0.1,50)
#plt.ylim(1e-5,0.9e4)
    
    # add label only for finest level data
    if ilev==n-1:
      plt.plot(x,np.log10(nh),"k-",label="$n_\mathrm{H}\, \left(\mathrm{cm}^{-3}\\right)$")
      plt.plot(x,np.log10(np.fabs(vr)),"b--",label="$\left|v_r\\right|\, \left(\mathrm{km\,s}^{-1}\\right)$")
      plt.plot(x,np.log10(tt),"r-.",label="Temperature (K)")
      plt.plot(x,np.log10(hp),"k:",label="$y(\mathrm{H}^{+})$")
    else:
      plt.plot(x,np.log10(nh),"k-")
      plt.plot(x,np.log10(np.fabs(vr)),"b--")
      plt.plot(x,np.log10(tt),"r-.")
      plt.plot(x,np.log10(hp),"k:")

  plt.legend(fontsize=12, loc="upper right")
  iy = str(i).zfill(5)
  opf = img_path+"/"+sim+"."+iy+".png"
  #plt.show()
  plt.savefig(opf, bbox_inches="tight")
  plt.close(fig)
  del fig
  dataio.close()
  del dataio
  del x, ro, nh, pg, vx, vr, tt, xh, xhe
  del xmin, xmax, dx, ng, time, D, c, n

quit()



