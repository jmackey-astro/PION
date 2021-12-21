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
    if(PION_NESTED_GRID)
        add_library(
            nested_grid

            setup_NG_grid.cpp
            setup_grid_NG_MPI.cpp
        )
    endif(PION_NESTED_GRID)
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
    if(PION_NESTED_GRID)
        add_library(
            nested_grid

            setup_NG_grid.cpp
        )
    endif(PION_NESTED_GRID)
endif(PION_PARALLEL)

target_link_libraries(grid PUBLIC boundaries)
target_link_libraries(grid PRIVATE global)
target_link_libraries(grid PRIVATE spdlog::spdlog)

if(PION_NESTED_GRID)
    target_link_libraries(nested_grid PUBLIC grid nested_boundaries)
    target_link_libraries(nested_grid PRIVATE global)
    target_link_libraries(nested_grid PRIVATE spdlog::spdlog)
endif(PION_NESTED_GRID)
