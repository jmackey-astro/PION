#!/bin/sh
#
# This is an example build script to compile PION on kay.ichec.ie
#
# (1)
# Make sure that you compiled the extra libraries with Gnu compilers...
# Check that line 33 of extra_libraries/install_all_libs.sh is *commented out*
# and reads '#KAY_INTEL=yes'
# Then run that script from the extra_libraries directory before building PION
# (2)
# Copy this file to the PION root directory and try building then.

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"
deps_dir="${HOME}/.local/"
build_dir="${script_dir}/build"

# load modules needed for PION
source /usr/share/Modules/init/bash
module load cmake3
module load openmpi/gcc

# Help Cmake find libraries if necessary
CMAKE_PREFIX_PATH="${deps_dir};${CMAKE_PREFIX_PATH}"
SUNDIALS_DIR="${deps_dir}"
SILO_DIR="${deps_dir}"

mkdir -p ${build_dir}
pushd ${build_dir}

cmake \
 -DCMAKE_BUILD_TYPE=Release      \
 -DPION_PARALLEL=ON              \
 -DPION_OMP=ON                   \
 -DPION_USE_SILO=ON              \
 -DPION_USE_FITS=OFF             \
 -DPION_UNIFORM_GRID=ON          \
 -DPION_NESTED_GRID=ON           \
 -DPION_TOOLS=ON                 \
 -DPION_BUILD_DOCUMENTATION=OFF  \
 -DPION_BUILD_TESTS=ON           \
 -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
 "${script_dir}"

make -j 20
popd

exit
