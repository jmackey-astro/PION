if(PION_PARALLEL)
    add_library(
        sim_control

        calc_timestep.cpp
        sim_control.cpp
        sim_control_MPI.cpp
        sim_init.cpp
        time_integrator.cpp
    )
    if(PION_NESTED_GRID)
    add_library(
        nested_sim_control

        sim_control_NG.cpp
        sim_control_NG_MPI.cpp
    )
    endif(PION_NESTED_GRID)
else()
    add_library(
        sim_control

        calc_timestep.cpp
        sim_control.cpp
        sim_init.cpp
        time_integrator.cpp
    )
    if(PION_NESTED_GRID)
        add_library(
            nested_sim_control

            sim_control_NG.cpp
        )
    endif(PION_NESTED_GRID)
endif(PION_PARALLEL)

if(PION_NESTED_GRID)
    target_link_libraries(nested_sim_control PUBLIC sim_control)
endif(PION_NESTED_GRID)

