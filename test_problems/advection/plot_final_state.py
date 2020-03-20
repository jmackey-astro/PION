#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import Silo
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/jm/code/pypion/Library")
import Plotting_Classes as ppion

file="advection_v020_t30_level00_0000.00007394.silo"
file="adv_l1_0000.00003697.silo"
plotting=ppion.Plotting2d(file)
n = plotting.nlevels()
plt.figure()

for ilev in range(n):
  if (n>1):
    s1 = "level"+ "%02d" % ilev
    if (file.find(s1) == -1):
      last=ilev-1
      s2 = "level"+ "%02d" % last
      file=file.replace(s2,s1)
      if (file.find(s1) == -1):
        print "error with filename, f=",file
        quit()

  plotting=ppion.Plotting2d(file)
  #d = plotting.reshaped_parameter2d("VelocityX")
  d = plotting.reshaped_parameter2d("Density")
  level_min = plotting.level_min()
  level_max = plotting.level_max()
  if (ilev==0):
    plt.xlim(level_min[0].value, level_max[0].value)
    plt.ylim(level_min[1].value, level_max[1].value)
  extents = [level_min[0].value, level_max[0].value, \
             level_min[1].value,level_max[1].value]
  plt.imshow(d,interpolation="nearest",cmap="viridis", \
             extent=extents,origin="lower", \
             vmin=0.0,vmax=10.0)
             #vmin=1.788839999,vmax=1.788840001)
  plotting.close()
  del plotting

plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.text(0.1,1.8,"Density, t=2")
#plt.show()
fn = "advection2D_l"+str(n)+".pdf"
plt.savefig(fn,bbox_inches="tight")
quit()
