cmake_minimum_required(VERSION 3.0.0)

project(
    pion
    VERSION 3.0.0
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
    PION_OMP
    "Build Pion with OpenMP. Default is disbled."
    OFF
)

option(
    PION_USE_SILO
    "Use the Silo library for data I/O. Default is enabled."
    ON
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

option(
    PION_TOOLS
    "Build postprocessing programmes in analysis dir.  Default is disabled."
    OFF
)

option(
    PION_SKIP_SOURCE
    "Do not build compile the Pion source files. Default is disabled."
    OFF
)

option(
    PION_BUILD_DOCUMENTATION
    "Build Pion documentation. Depends on Doxygen."
    OFF
)

option(
    PION_INTEL
    "Build with Intel compilers and mathimf.h math library"

    OFF
)

#  This is required on OSX so that the Boost code will compile successfully
set(CMAKE_CXX_STANDARD 14)

if(NOT PION_SKIP_SOURCE)
    if(PION_PARALLEL)
        add_definitions(-DPARALLEL)
        add_definitions(-DUSE_MPI)
    else()
        add_definitions(-DSERIAL)
    endif(PION_PARALLEL)

    if(PION_USE_SILO AND PION_USE_FITS)
        message("WARNING: building with both Silo and Fits libraries compiled...")
    endif(PION_USE_SILO AND PION_USE_FITS)

    if(NOT PION_USE_SILO AND NOT PION_USE_FITS)
        message("Cannot build without Silo or Fits libraries...")
        message("Please rerun the CMake command with one of them set to `ON`.")
        message("Eg. -DPION_USE_FITS=ON.")
        return()
    endif(NOT PION_USE_SILO AND NOT PION_USE_FITS)

    if(PION_USE_SILO)
        add_definitions(-DSILO)
    endif(PION_USE_SILO)

    if(PION_USE_FITS)
        add_definitions(-DFITS)
    endif(PION_USE_FITS)

    if(PION_OMP)
        add_definitions(-DPION_OMP)
    endif(PION_OMP)

    if(PION_INTEL)
        add_definitions(-DINTEL)
    endif(PION_INTEL)


    #
    # Find dependant libraries
    #
    set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

    find_package(Boost 1.73.0 REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})

    find_package(SUNDIALS REQUIRED)

    add_definitions("-DCVODE${SUNDIALS_VERSION_MAJOR}")

    if(PION_USE_SILO)
        find_package(SILO REQUIRED)
        include_directories(${SILO_INCLUDE_DIR})
    endif(PION_USE_SILO)

    if(PION_USE_FITS)
        find_package(FITS REQUIRED)
        find_package(ZLIB REQUIRED)
        find_package(CURL REQUIRED)
        find_package(BZip2 REQUIRED)
        find_package(Threads REQUIRED)
        include_directories(${FITS_INCLUDE_DIR})
    endif(PION_USE_FITS)

    if(PION_PARALLEL)
        find_package(MPI REQUIRED COMPONENTS CXX)
        include_directories(${MPI_CXX_INCLUDE_DIRS})
    endif(PION_PARALLEL)

    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        find_package(Curses)
        find_package(READLINE REQUIRED)
        include_directories(${READLINE_INCLUDE_DIR})
    endif(CMAKE_BUILD_TYPE MATCHES "Debug")

    if (PION_OMP)
        find_package(OpenMP REQUIRED COMPONENTS CXX)
        include_directories(${OpenMP_CXX_INCLUDE_DIR})
    endif(PION_OMP)

    find_package(fmt REQUIRED)
    find_package(spdlog REQUIRED)

    add_definitions("-DSPDLOG_FMT_EXTERNAL")
    add_definitions("-DSPDLOG_COMPILED_LIB")


    if (${spdlog_VERSION_MAJOR} GREATER_EQUAL 1)
        if (${spdlog_VERSION_MINOR} GREATER_EQUAL 6)
            message(STATUS "spdlog version = ${spdlog_VERSION_MAJOR}.${spdlog_VERSION_MINOR}")
            add_definitions("-DSPDLOG_FWD")
        endif (${spdlog_VERSION_MINOR} GREATER_EQUAL 6)
    endif (${spdlog_VERSION_MAJOR} GREATER_EQUAL 1)

    #
    # Include optional subdirectories
    #

    include_directories(
        SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source
    )
endif(NOT PION_SKIP_SOURCE)


#
# Setup the Pion library
#
if(NOT PION_SKIP_SOURCE)
    add_subdirectory(source)
endif(NOT PION_SKIP_SOURCE)

if(PION_BUILD_DOCUMENTATION)
    add_subdirectory(doc)
endif(PION_BUILD_DOCUMENTATION)

if(PION_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
    set(PION_TOOLS ON)
endif(PION_BUILD_TESTS)

if(PION_TOOLS)
    add_subdirectory(analysis)
    if(NOT PION_PARALLEL)
      message("WARNING: Requested PION_TOOLS: most tools need PION_PARALLEL")
      message("WARNING: but this is not set.  Most tools will not be compiled")
    endif(NOT PION_PARALLEL)
endif(PION_TOOLS)

# Output project build options
message(STATUS "Build options:")
message(STATUS "    CMAKE_BUILD_TYPE         = ${CMAKE_BUILD_TYPE}")
message(STATUS "    PION_PARALLEL            = ${PION_PARALLEL}")
message(STATUS "    PION_OMP                 = ${PION_OMP}")
message(STATUS "    PION_INTEL               = ${PION_INTEL}")
message(STATUS "    PION_USE_SILO            = ${PION_USE_SILO}")
message(STATUS "    PION_USE_FITS            = ${PION_USE_FITS}")
message(STATUS "    PION_NESTED_GRID         = ${PION_NESTED_GRID}")
message(STATUS "    PION_UNIFORM_GRID        = ${PION_UNIFORM_GRID}")
message(STATUS "    PION_TOOLS               = ${PION_TOOLS}")
message(STATUS "    PION_SKIP_SOURCE         = ${PION_SKIP_SOURCE}")
message(STATUS "    PION_BUILD_TESTS         = ${PION_BUILD_TESTS}")
message(STATUS "    PION_BUILD_DOCUMENTATION = ${PION_BUILD_DOCUMENTATION}")
