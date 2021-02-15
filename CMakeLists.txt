cmake_minimum_required(VERSION 3.0.0)

project(
    pion
    VERSION 0.0.1
    LANGUAGES CXX
)

#
# Set subdirectory options
#


#
# Setup installation path options
#

# Sane default install directories
include(GNUInstallDirs)

# User provided install directories
set(
    PION_INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE STRING
    "The install directory for Pion executables."
)

set(
    PION_INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE STRING
    "The install directory for Pion libraries."
)

set(
    PION_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE STRING
    "The install directory for Pion headers."
)

set(
    PION_INSTALL_DATADIR ${CMAKE_INSTALL_DATADIR} CACHE STRING
    "The install directory for Pion data (eg. cmake scripts)."
)


#
# Set build type options
#

option(
    PION_PARALLEL
    "Build Pion with MPI. Default is enabled."
    ON
)

option(
    PION_USE_SILO
    "Use the Silo library for data I/O. Default is enabled."
    OFF
)

option(
    PION_USE_FITS
    "Use the FITS library for data I/O. Default is disabled."
    OFF
)

option(
    PION_NESTED_GRID
    "Build the nested grid simulation(s). Default is enabled."
    ON
)

option(
    PION_UNIFORM_GRID
    "Build the uniform grid simulation(s). Default is disabled."
    OFF
)

if(PION_PARALLEL)
    add_definitions(-DPARALLEL)
    add_definitions(-DUSE_MPI)
else()
    add_definitions(-DSERIAL)
endif(PION_PARALLEL)

if(PION_USE_SILO AND PION_USE_FITS)
    message("Cannot build with both Silo and Fits libraries...")
    message("Please rerun the CMake command without one of them set to `ON`.")
    return()
endif(PION_USE_SILO AND PION_USE_FITS)

if(PION_USE_SILO)
    add_definitions(-DSILO)
endif(PION_USE_SILO)

if(PION_USE_FITS)
    add_definitions(-DFITS)
endif(PION_USE_FITS)

#
# Find dependant libraries
#

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(GSL REQUIRED COMPONENTS spline)
find_package(SUNDIALS REQUIRED)

add_definitions("-DCVODE${SUNDIALS_VERSION_MAJOR}")

if(PION_USE_SILO)
    find_package(SILO REQUIRED)
endif(PION_USE_SILO)

if(PION_USE_FITS)
    find_package(FITS REQUIRED)
    find_package(ZLIB REQUIRED)
    find_package(CURL REQUIRED)
    find_package(BZip2 REQUIRED)
    find_package(Threads REQUIRED)
endif(PION_USE_FITS)

if(PION_PARALLEL)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif(PION_PARALLEL)


#
# Include optional subdirectories
#

include_directories(
    SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source
)

include_directories(
    SYSTEM PUBLIC ${SILO_INCLUDE_DIR}
)

#
# Setup the Pion library
#

add_subdirectory(source)

# Output project build options
message(STATUS "Build options:")
message(STATUS "    PION_PARALLEL     = ${PION_PARALLEL}")
message(STATUS "    PION_USE_SILO     = ${PION_USE_SILO}")
message(STATUS "    PION_USE_FITS     = ${PION_USE_FITS}")
message(STATUS "    PION_NESTED_GRID  = ${PION_NESTED_GRID}")
message(STATUS "    PION_UNIFORM_GRID = ${PION_UNIFORM_GRID}")
