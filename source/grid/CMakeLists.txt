add_library(
    grid

    cell_interface.cpp
    setup_fixed_grid.cpp
    setup_fixed_grid_MPI.cpp
    setup_grid_NG_MPI.cpp
    setup_NG_grid.cpp
    stellar_wind_angle.cpp
    stellar_wind_BC.cpp
    uniform_grid.cpp
    uniform_grid_pllel.cpp
)

target_link_libraries(grid PUBLIC 
    dataio
    microphysics
    raytracing
    spatial_solvers
)
