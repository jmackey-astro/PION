#!/bin/sh

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"
deps_dir="${script_dir}/extra_libraries"
build_dir="${script_dir}/build"

mkdir -p ${build_dir}
pushd ${build_dir}

CMAKE_PREFIX_PATH="${deps_dir};${CMAKE_PREFIX_PATH}"
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DPION_USE_SILO=ON \
    -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
    -DPION_UNIFORM_GRID=ON \
    -DPION_NESTED_GRID=ON \
    -DPION_TOOLS=ON \
    -DPION_OMP=OFF \
    -DPION_BUILD_TESTS=ON \
    "${script_dir}"

make -j 4
popd
