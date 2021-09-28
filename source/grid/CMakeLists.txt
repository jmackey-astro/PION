if(PION_PARALLEL)
    add_library(
        grid

        cell_interface.cpp
        setup_fixed_grid.cpp
        setup_fixed_grid_MPI.cpp
        stellar_wind_angle.cpp
        stellar_wind_BC.cpp
        stellar_wind_latdep.cpp
        uniform_grid.cpp
        uniform_grid_pllel.cpp
    )
    add_library(
        nested_grid

        setup_NG_grid.cpp
        setup_grid_NG_MPI.cpp
    )
else()
    add_library(
        grid

        cell_interface.cpp
        setup_fixed_grid.cpp
        stellar_wind_angle.cpp
        stellar_wind_BC.cpp
        stellar_wind_latdep.cpp
        uniform_grid.cpp
    )
    add_library(
        nested_grid

        setup_NG_grid.cpp
    )
endif(PION_PARALLEL)


target_link_libraries(grid PUBLIC 
    dataio
    microphysics
    raytracing
    spatial_solvers
)

target_link_libraries(nested_grid PUBLIC grid)

if(PION_OMP)
    target_link_libraries(grid PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(nested_grid PUBLIC OpenMP::OpenMP_CXX)
endif(PION_OMP)

