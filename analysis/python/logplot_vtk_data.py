# -*- coding: iso-8859-15 -*-

# Author: Jonathan Mackey
# Date:   2016.06.[20-22]
#
# Description:
# - For a list of VTK files, the sript will
# - read in VTK data produced by PION 2D projection code.
# - plot the column density and nebular line emission on log scale.
#
# Needs python-vtk, python-astropy, python-numpy, python-matplotlib

#####################################################################
# plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

# numpy
import numpy as np

# vtk
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN

# parse arguments
import os
import argparse

# convert from CGS units to others:
import astropy.constants as apc
from astropy import units as u

# Use LaTeX in figure labels:
from matplotlib import rc
rc('text', usetex=True)
#####################################################################

#####################################################################
# Parse arguments and get list of files to read
#
parser = argparse.ArgumentParser(description='Process some values.',
 usage='logplot_vtk_data.py <path-to-files> <base-filename> <image-path> <image-filename>')
parser.add_argument('Path', help='path to VTK files')
parser.add_argument('file_base', help='base filename of VTK files')
parser.add_argument('img_path', help='path to save images into')
parser.add_argument('img_base', help='base filename of image files')
args = parser.parse_args()
#print args

# make sure that path has a trailing /
file_path = args.Path
if not file_path.endswith('/'):
  file_path += "/"

img_path  = args.img_path
img_file  = args.img_base

#Reads all the files in a directory and then puts their names into a 1D array.
data_files = [f for f in os.listdir(file_path)] # if isfile(join(file_path, f))]
# Sort list by timestep
data_files.sort()

# Remove non-VTK files, and non-requested files.
search = ".vtk"
data_files = [f for f in data_files if search in f]
search = args.file_base
data_files = [f for f in data_files if search in f]

# Add path to files
data_files_with_path = [f.replace(f,file_path+f) for f in data_files]
#####################################################################

VTKFiles = len(data_files_with_path)
print "%s VTK files detected." %VTKFiles
Figure = 0

#####################################################################
# Loop over all files founds in directory matching the search string.
for filename in data_files_with_path:
  print "Reading from file: ",filename
  reader = vtkStructuredPointsReader()
  reader.SetFileName(filename)
  reader.ReadAllScalarsOn()
  reader.Update()
  data = reader.GetOutput()

# This gives the number of vertices in each dimension
  dim = data.GetDimensions()

# This gives the lower and upper bounds of vertex number in each dimension
  bounds = data.GetBounds()
  bounds = bounds[0:4]
  bounds = np.array(bounds)
  bounds = bounds/apc.pc.cgs.value
# bounds1 is for setting extents of the lower half-plane (reflected).
  bounds1 = [bounds[0],bounds[1],-bounds[3],-bounds[2]]

# This gives the extents of simulation (Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
  extent = data.GetExtent()

# Get time and cycle (PION-SPECIFIC CODE)
  fdata = data.GetFieldData()
  time = fdata.GetAbstractArray(0).GetValue(0) *u.second
  time = time.to(u.Myr)
  cycle = fdata.GetAbstractArray(1).GetValue(1)

  # Sigma is the projected mass per cm^2
  Sigma  = VN.vtk_to_numpy(data.GetCellData().GetArray('AA_SurfaceMass'))
  # Halpha and NII are the respective line intensities at each pixel,
  # in erg/cm^2/s/sq.arcsec.
  Halpha = VN.vtk_to_numpy(data.GetCellData().GetArray('AA_Halpha'))
  NII    = VN.vtk_to_numpy(data.GetCellData().GetArray('AA_NII_ll6584'))

  dim = np.array([extent[1],extent[3]])

# Reshape and edit Sigma and Lines arrays for plotting.
# 'F' means Fortran ordering of array (this is what works...).
  Sigma = Sigma.reshape(dim, order='F')
  Sigma = Sigma.T
  Sigma = np.flipud(Sigma)
  Sigma /= 2.338e-24  # Convert to H column density for ISM mass frac
  Sigma = np.log10(Sigma)

  Lines = Halpha+NII
  Lines = Lines.reshape(dim, order='F')
  Lines = Lines.T
  Lines *= 1.767e17   # Convert to Rayleigh (at Halpha wavelength)
  Lines = np.log10(Lines)


# Set up figure.
  fig = plt.figure()
  ax1 = fig.add_subplot(2,1,1)
  ax1.set_title("Time = %5.4g Myr"%time.value)

# Add column density plot to upper half-plane
  im1 = ax1.imshow(Sigma,interpolation = 'nearest', cmap = 'Greys', extent=bounds)
  plt.colorbar(im1, format="%4.2f")
  ax1.axes.get_xaxis().set_visible(False)
  ax1.axes.set_ylabel("R (pc)")
  ax1.text(0.75,0.92,'N(H) (cm$^{-2}$)', transform=ax1.transAxes)

# Add nebular-line-emission plot to lower half-plane.
  ax2 = fig.add_subplot(2,1,2)
  im2 = ax2.imshow(Lines,interpolation = 'nearest', cmap = 'Reds', extent=bounds1)
  norm = mpl.colors.Normalize(0, vmax=Lines.max())
  im2.set_norm(norm) 
  plt.colorbar(im2, format="%4.2f")
  ax2.text(0.75,0.05,'H$\\alpha$ (Rayleigh)', transform=ax2.transAxes)
  ax2.axes.set_xlabel("Z (pc)")

# Save figures
  fig.subplots_adjust(wspace=0, hspace=0)
  nn = "%03d" % (Figure)
  outfile="%s/%s_%s.png" %(img_path,img_file,nn)
  print "   Saving to file:",outfile
  plt.savefig(outfile, bbox_inches='tight', dpi=300)
  plt.close(fig)
# increment counter and continue loop.
  Figure += 1
# End of loop over all files
#####################################################################

exit()

