# Makefile for my Grid Hydro Code, Jonathan Mackey
SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cc .o
srcdir = .

CC = g++

CXXFLAGS = -Wall #-O2 -g
# --------------------
# Compile-time options
# --------------------
# --dimensionality of grid; only one can be uncommented!!!
#OPT += -DGRID_1D
OPT += -DGRID_2D
#OPT += -DGRID_3D

# Type of Grid
OPT += -DUNIFORM_FV

# Number of space dimensions for equations
# NOT USING THIS ANYMORE -- HARD-CODED SET TO 3 FOR SIMPLICITY.
#OPT += -DEQN_1D
#OPT += -DEQN_2D
#OPT += -DEQN_3D

# What equations/variables to solve
OPT += -DHD
#OPT += -DMHD

# Type of Solver
OPT += -DLAX_FRIEDRICHS
#OPT += -DGODUNOV

# Type of Initial Condition
#OPT += -DSHOCKTUBE

CXXFLAGS += $(OPT)

H5include = -I/usr/include 
H5libs =  -lz -lm -lpthread /usr/lib/libhdf5_cpp.a /usr/lib/libhdf5.a
FITSlibs = -lcfitsio
#H5flags = -fomit-frame-pointer

h5tofits: h5tofits.o
	$(CC) $(CXXFLAGS) $(H5include) h5tofits.o -o h5tofits $(H5libs) $(FITSlibs)

h5tofits.o: h5tofits.cc

clean:
	rm -f *.o h5tofits
