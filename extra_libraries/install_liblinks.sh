#!/bin/bash

# If you already have pion libraries installed somewhere else, this
# file will link to them, just edit the path.

#mkdir include
#mkdir bin
#mkdir lib

#./install_cfitsio.sh
#./install_cvodes.sh
#./install_silo.sh

lib_path=../../pion/extra_libraries
ln -s ${lib_path}/lib ./lib
ln -s ${lib_path}/include ./include
ln -s ${lib_path}/bin ./bin

