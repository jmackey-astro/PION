add_library(
    spatial_solvers

    solver_eqn_base.cpp
    solver_eqn_hydro_adi.cpp
    solver_eqn_mhd_adi.cpp
)

target_link_libraries(spatial_solvers PUBLIC 
    riemann_solvers
    coord_sys
)
target_link_libraries(spatial_solvers PRIVATE spdlog::spdlog)
