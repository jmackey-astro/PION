#!/bin/sh

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"
build_dir="${script_dir}/build"

mkdir -p ${build_dir}
pushd ${build_dir}

cmake \
 -DCMAKE_BUILD_TYPE=Release      \
 -DCMAKE_CXX_COMPILER=mpicxx     \
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

make -j 4
popd 

