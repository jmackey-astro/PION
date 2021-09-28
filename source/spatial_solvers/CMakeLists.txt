add_library(
    spatial_solvers

    solver_eqn_base.cpp
    solver_eqn_hydro_adi.cpp
    solver_eqn_mhd_adi.cpp
)

target_link_libraries(spatial_solvers PUBLIC 
    riemann_solvers
    equations
    coord_sys
)

if(PION_OMP)
    target_link_libraries(spatial_solvers PUBLIC OpenMP::OpenMP_CXX)
endif(PION_OMP)
