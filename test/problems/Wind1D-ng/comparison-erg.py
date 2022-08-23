# -*- coding: iso-8859-15 -*-
''' This code compares the mass of the stellar wind in MESA and in the PION run'''
    
#pion imports:
#base_path='/vol/aibn128/data1/yfichtner/'
#base_path='/mnt/massive-stars/share/pypion/'
base_path='/home/jm/.local'
import sys
sys.path.insert(0,base_path+"/silo/lib")
import Silo
from pypion.ReadData import ReadData
import glob
from astropy import units as u
from astropy import constants as apc

#import mesa_reader as mr

#general imports
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('classic')
import matplotlib as mpl
mpl.rcParams["font.size"]=20
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['lines.linewidth'] = 1.5
mpl.use('Agg')
#plt.use('Agg')
import os
os.environ["OMP_NUM_THREADS"] = "4"
from sys import version_info
version_info >= (3, 7)

################
#path="/vol/aibn128/data1/yfichtner/Pion_mesa_new/Test_Version2"
#subdirectory="Test_1_GALrotation_1.700_0.0_standard"
path="silo"
energies = []
times = []
for sim in ["Galrot_test_1.700_0.0_d1l2n0128", "Galrot_test_1.700_0.0_d1l1n0256", "Galrot_test_1.700_0.0_d1l1n0512", "Galrot_test_1.700_0.0_d1l3n0128", "Galrot_test_1.700_0.0_d1l4n0128", "Galrot_test_1.700_0.0_d1l5n0128"]:

##############

#### Get pion values

  #files=[]
  #for file_h in glob.glob("%s/Galrot_test_1.700_0.0_level.silo" %(path,sim)):
  #        files.append(file_h)	#maybe need sorted
  #files = sorted(files)
  if sim.find("d1l2") > 0:
    files0 = sorted(glob.glob(path+"/"+sim+"_level00_0000.*.silo"))
    files1 = sorted(glob.glob(path+"/"+sim+"_level01_0000.*.silo"))
    lev=2
  elif sim.find("d1l1") > 0:
    files0 = sorted(glob.glob(path+"/"+sim+"_0000.*.silo"))
    lev=1
    print(files0)
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
    print("can't find requested files")
    quit()

  erg_p = []		#total thermal and kinetic energy in the box
  t = []

  for i in np.arange(0,len(files0),1):
    #print(i)
    if lev==1:
      print(files0[i])
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
    
    D = dataio.get_1Darray("Density")
    density  = D['data'] # g/c^3
    vr       = (dataio.get_1Darray("VelocityX")['data'])  #cm/s
    erg      = (dataio.get_1Darray("InternalEnergy")['data']) # erg/g
    if lev>1:
      vM       = dataio.get_1Darray("NG_Mask")['data']
    sim_time = (dataio.sim_time()).to(u.Myr)
    xmax     = D['max_extents']
    xmin     = D['min_extents']
    ng       = dataio.ngrid()
    n        = dataio.nlevels()
    dataio.close()
    
    energy_tot = 0.0
    for ilev in range(n):
      print(ilev,files,xmin,xmax)
      lim_min=xmin[ilev][0]
      lim_max=xmax[ilev][0]
      dx=(lim_max-lim_min)/ng[0]
      rmin = np.arange(lim_min   ,lim_max,    dx)
      rmax = np.arange(lim_min+dx,lim_max+dx, dx)
      m      = 4.0*np.pi/3.0 * (rmax**3 - rmin**3) # vol
      if lev>1:
        energy = m * vM[ilev]* (0.5*density[ilev]*vr[ilev]**2 + density[ilev]*erg[ilev])
      else:
        energy = m * (0.5*density[ilev]*vr[ilev]**2 + density[ilev]*erg[ilev])
      eee    = np.sum(energy)
      energy_tot += eee

    t.append(sim_time.value)
    erg_p.append(energy_tot)

  energies.append(erg_p)
  times.append(t)

print(len(t),len(energies[0]),len(energies[1]))

# theoretical input power of wind
power = 0.5 * (1.0e-6*(apc.M_sun.cgs).value*1e6) * 2.5e8**2

theory=[]
for j in range(len(energies)):
  energies[j] = np.array(energies[j])
  times[j] = np.array(times[j])
  theory.append(energies[j][0]+power*np.array(times[j]))

## MASS comparison
plt.figure()
plt.plot(np.array(times[0]),(energies[0]-theory[0])/theory[0],":",label='l2 n128')
plt.plot(np.array(times[1]),(energies[1]-theory[1])/theory[1],"--",label='l1 n256')
plt.plot(np.array(times[2]),(energies[2]-theory[2])/theory[2],"--",label='l1 n512')
plt.plot(np.array(times[3]),(energies[3]-theory[3])/theory[3],":",label='l3 n128')
plt.plot(np.array(times[4]),(energies[4]-theory[4])/theory[4],":",label='l4 n128')
plt.plot(np.array(times[5]),(energies[5]-theory[5])/theory[5],":x",label='l5 n128')
#plt.plot(np.array(times[2]),np.array(energies[2]),"-",label='Pion: Total energy 512')
#plt.plot(np.array(times[3]),np.array(energies[3]),"-.",label='Pion: Total energy 1024')
#plt.plot(np.array(times[4]),np.array(energies[4]),"-",label='Pion: Total energy 2048')

plt.xlabel('Time [Myr]')
#plt.xlim(0,0.5)
plt.ylabel('relative error [total energy]')
#plt.ylim(0,0.5)
plt.legend(loc='upper left',fontsize=14,ncol=2,frameon=False,labelspacing=0.25,columnspacing=1)
plt.tight_layout()
plt.savefig('Comparison_erg.png')
#  plt.show()


quit()
