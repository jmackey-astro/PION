# PION

PION is open source software primarily for simulation of nebulae around stars, but also for other astrophysical applications.  It has an open development model, where anyone can download and use the software, change it to suit their needs, and contribute changes back to the user community.


## Information

 * The PION homepage is at [https://www.pion.ie](https://www.pion.ie), with documentation including a quick-start guide at [https://www.pion.ie/docs/](https://www.pion.ie/docs/).

 * The git repository to obtain the latest version of PION is at [https://git.dias.ie/massive-stars-software/pion](https://git.dias.ie/massive-stars-software/pion)

 * This is mirrored on github at [https://github.com/jmackey-astro/PION](https://github.com/jmackey-astro/PION)

 * Python routines and libraries for visualising results are at [https://git.dias.ie/massive-stars-software/pypion](https://git.dias.ie/massive-stars-software/pypion)

 * Contact [info@pion.ie](mailto:info@pion.ie) for help and information


## Build Instructions

To configure and build Pion, a C++14 compatible compiler is required, and CMake >= 3.0. For developers, it is recommended to create a script similar to the following example in the root directory of the project for repeated building of the project. 

```
#!/bin/sh

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"
#deps_dir="${script_dir}/extra_libraries"
#modules_path="/path/to/modules"
build_dir="${script_dir}/build"

# (optionally) load any modules
#module load cmake3
#module load gcc
#module load gsl/gcc
#module load openmpi/gcc/4.0.5

# Help Cmake find libraries if necessary
#CMAKE_PREFIX_PATH="${modules_path}/gsl/gcc/2.5;${deps_dir};${CMAKE_PREFIX_PATH}"
#SUNDIALS_DIR="${deps_dir}/share/sundials"
#SILO_DIR="${deps_dir}"

mkdir -p ${build_dir}
pushd ${build_dir}

cmake \
    #-DCMAKE_BUILD_TYPE=Debug \
    #-DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
    #-DSUNDIALS_DIR="${SUNDIALS_DIR}" \
    #-DUSE_SILO=ON \
    #-DSILO_DIR="${SILO_DIR}" \
    [more options below] \
    "${script_dir}"

make -j 4

popd # build_dir
```

An example script for [Debian10/Ubuntu20.04 is here](https://homepages.dias.ie/jmackey/pion-dev-doc/_downloads/c09c3b85afae16c99963dfd1d25965d5/build_debian.sh).
To install execute `make -C ${build_dir} install`.

### Useful Options

#### Setting Compile/Link Flags
Option | Effect
------ | ------
`-DCMAKE_BUILD_TYPE=...`   | Set to `Release` or `Debug` for release or debug builds. Sets a number of flags by default.
`-DCMAKE_CXX_COMPILER=...` | Set the C++ compiler.
`-DCMAKE_CXX_FLAGS=...`    | Set the flags to pass to the C++ compiler explicitly. Overrides the default flags.

#### Enabling/Disabling Aspects of Pion
Option | Effect
------ | ------
`-DPION_NESTED_GRID=...`          | Set to `ON` to build nested grid simulations (default `ON`).
`-DPION_UNIFORM_GRID=...`         | Set to `ON` to build uniform grid simulations (default `OFF`).
`-DPION_PARALLEL_=...`            | Set to `ON` to enable MPI for the Pion build, or `OFF` to disable MPI (default `ON`).
`-DPION_USE_SILO=...`             | Set to `ON` to use Silo for handling data I/O (default `OFF`).
`-DPION_USE_FITS=...`             | Set to `ON` to use Fits for handling data I/O (default `OFF`).
`-DPION_TOOLS=...`                | Set to `ON` to also compile support programs in analysis subdir (default `OFF`).
`-DPION_SKIP_SOURCE=...`          | Set to `ON` to skip compilation of the Pion source files (default `OFF`).
`-DPION_BUILD_DOCUMENTATION =...` | Set to `ON` to build Doxygen documentation (default `OFF`).
`-DPION_BUILD_TESTS=...`          | Set to `ON` to build tests (default `OFF`).


#### Dependency Helpers
Option | Effect
------ | ------
`-DSUNDIALS_DIR=...` | Instruct CMake to search the provided path for the Sundials library.
`-DSILO_DIR=...` | Instruct CMake to search the provided path for the Silo library.
`-DFITS_DIR=...` | Instruct CMake to search the provided path for the Fits library.
`-DBOOST_ROOT=...` | Instruct CMake to search the provided path for the Boost library.

#### Setting installation directories
Option | Effect
------ | ------
`-DCMAKE_INSTALL_PREFIX=...` | Set the root install directory for the compiled libraries and programs. 
`-DPION_INSTALL_BINDIR=...`  |  Set the install directory for Pion executables. Use a relative path to set the path relative to `${CMAKE_INSTALL_PREFIX}` (default `bin`).
`-DPION_INSTALL_LIBDIR=...`  |  Set the install diretory for Pion libraries. Use a relative path to set the path relative to `${CMAKE_INSTALL_PREFIX}` (default `lib`).
`-DPION_INSTALL_INCLUDEDIR=...`  |  Set the install diretory for Pion header files. Use a relative path to set the path relative to `${CMAKE_INSTALL_PREFIX}` (default `include`).
`-DPION_INSTALL_DATADIR=...`  |  Set the install diretory for Pion data (eg. CMake scripts). Use a relative path to set the path relative to `${CMAKE_INSTALL_PREFIX}` (default `share`).

### Documentation 

The Doxygen documentation can be built with `make -C [build-dir] doc` provided the `-DPION_BUILD_DOCUMENTATION=ON` flag was provided to CMake.

## Developers

The following people have contributed to the development of PION:

  * Harpreet Dhanoa
  * Margueritta Goulden   
  * Samuel Green
  * Robert Kavanagh
  * Andrew Lim
  * Jonathan Mackey       <jmackey@cp.dias.ie>
  * Maria Moutzouri    
  * Ciarán O'Rourke
  * Davit Zargaryan


## License

PION is distributed under a BSD3 License.  Downloading, using, modifying and/or re-distributing the software implies acceptance of the License.  See [LICENSE.md](https://git.dias.ie/massive-stars-software/pion/-/blob/master/LICENSE.md).
