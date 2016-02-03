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
#H5flags = -fomit-frame-pointer

twod2oned45deg: twod2oned45deg.o
	$(CC) $(CXXFLAGS) $(H5include) twod2oned45deg.o -o twod2oned45deg $(H5libs)

twod2oned45deg.o: twod2oned45deg.cc

clean:
	rm -f *.o main
