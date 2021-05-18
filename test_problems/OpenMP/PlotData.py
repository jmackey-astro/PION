# Author: Sam Green, Created: 11-10-20

# Script to call all PyPion classes to plot Silo data.

# JM: added paths to PyPion and Python-Silo
# make sure you have activated modules:
# $ module load intel/2018u4
# $ module load conda
# $ source activate
# 
root_path="/ichec/home/users/jmackey/active"
import sys
sys.path.append(root_path+"/pion/extra_libraries/python/lib")
import Silo
sys.path.append(root_path+"/pion_python/Library")

from Plotting_Classes import Plotting2d, Plotting3d
from argparse_command import InputValues
#-------------------------------
import matplotlib
# Using this to stop matplotlib from using a $DISPLAY environment variable.
# i.e. This now works over ssh without the need for a x-server.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#-------------------------------
import warnings
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)
#-------------------------------
line = InputValues()
time_dicts = line.time_dicts
dimen = line.dimen
#-------------------------------

for files in time_dicts:

	arr = time_dicts[files]

	var1 = ["Density", -22, -27, "viridis", 'y', 63]
	# var1 = ["Temperature", 8, 3, "inferno", 'log', 'y', 127]

	fig = plt.figure()

        if dimen == "2d" or dimen == "2D":

             # print(dimen)
             a = Plotting2d(arr).plot2d_1(var1[0], fig, var1)

        elif dimen == "3d" or dimen == "3D":

             # print(dimen)
             a = Plotting3d(arr).XZXYslice(var1[0], fig, var1)

        else:
             print("Please choose a correct dimenion (1d, 2d, or 3d)")

	imagefile = "%s%s_%s.png" % (line.img_path, line.img_file, time_dicts[files][0][len(time_dicts[files][0]) - 13:len(time_dicts[files][0]) - 5])
	plt.savefig(imagefile, bbox_inches='tight', dpi=300)

