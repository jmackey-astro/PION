# Sam Green 2016
# This script reads in .silo files and creates contour plots of the data inside.

# -2016-05-24 SG: Created to work with the run_double_Mach_reflection_test.sh script in the pion code.



import os
import Silo
import numpy as np
import matplotlib.pyplot as plt
import pylab

from os import listdir
from os.path import isfile, join

if os.access("multi_ucd3d.pdb",os.R_OK):
    file = "multi_ucd3d.pdb"
elif os.access("multi_ucd3d.h5",os.R_OK):
    file = "multi_ucd3d.h5"


#Reads all the files in a directory and then puts their names into a 1D array.
file_path = "SILO/"
data_files = [f for f in listdir(file_path) if isfile(join(file_path, f))]

print "The figures are now being made!"

for files in data_files:
	#Opens and reads the .silo files.
	data = "SILO/%s" %files
	db = Silo.Open(data)


#	print "db = ",db
#	print "db.filename = '%s'"%db.filename

	toc = db.GetToc()
	print "\n-- TOC --\n",toc

#	print "cycle='%d'"%db.GetVar("cycle")
#	print "dtime='%f'"%db.GetVar("dtime")
#	print "_fileinfo='%s'"%db.GetVar("_fileinfo")
	
	#Reads in the dim array and then flips the dimensions around.
	dims = db.GetVar("Density_dims")
	dims = dims[::-1]
	
	#Reads in the relevant arrays.
	density = db.GetVar("Density_data")
	max_extents = db.GetVar("UniformGrid_max_extents")
	min_extents = db.GetVar("UniformGrid_min_extents")
	
	#Define the max and min values of the x and y axis.
	xmax = max_extents[0]
	xmin = min_extents[0]
	ymax = max_extents[1]
	ymin = min_extents[1]
	
	#Puts the array into the correct format, i.e. (a,b)
	density = np.array(density).reshape(dims)
	
	#Creates the the string for the plot titles
	Title = "%s" %files
	if Title[19:23] == "Hcor":
		title = "Nx=%s, %s, %s, Timestep=%s"%(Title[11:14], Title[15:18], Title[19:23], Title[29:32])
	else:
		title = "Nx=%s, %s, %s, Timestep=%s"%(Title[11:14], Title[15:18], Title[19:21], Title[29:32])	

	#Plots the arrays as contour plots. The extent function scales the x and y axis to a specific set of values
	plt.figure(figsize=(10,3.5))
	plt.contour(density, extent=[xmin, xmax, ymin, ymax])
	plt.text(1, 0.9, '%s'%title)


	plt.savefig("FIGS/%s.png" %files)
