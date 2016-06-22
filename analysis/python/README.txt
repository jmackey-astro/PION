
---------------------------------------------------------------------
This is a directory for python scripts that are useful for postprocessing PION simulations.  Most of them are hardcoded to tackle specific problems, so you will probably have to edit them to make them work for you.
---------------------------------------------------------------------

---------------------------------------------------------------------
plot_vtk_data.py:
logplot_vtk_data.py:

These files read in VTK datafiles produced by
PION/analysis/Project2D_data/project_data2D  and plot neutral H column density
and H-alpha+[NII] emission on the upper and lower half-planes, respectively.
The mean mass per H atom is hard-coded, there is assumed to be no molecular
phase, and it is assumed that the ionization state of N is the same as that of
H.  Linear scales are used for plot_vtk_data.py, whereas logarithmic scales are
used for logplot_vtk_data.py.
---------------------------------------------------------------------





