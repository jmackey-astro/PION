Directories:

comms: the parallel communication classes.

coord_sys: The VectorOps classes, for calculating coordinate-dependent vector operations like div, curl, grad, and the centre-of-volume of the cell.

equations: The systems of equations to solve.

Riemann_solvers: functions to calculate inter-cell fluxes across boundaries, including some solvers that are not Riemann solvers!

flux_calc: a class a level above the Riemann solvers, controlling which solver is called, and including the fluxes of the tracer variables (not included in Riemann solvers).

spatial_solvers: this is for the grid-aware classes which calculate fluxes across the grid.  Maybe they should be in the grid/ subdir?  but this is ok.

grid: data structures and access functions for the cell and uniform_grid classes.

dataIO: classes for data input and output.

defines: include files with some #defined variables controlling the simulation.

ics: initial condition generation classes.

microphysics: microphysics classes.

raytracing: raytracing classes

