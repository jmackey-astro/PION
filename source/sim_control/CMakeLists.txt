if(PION_PARALLEL)
    add_library(
        sim_control

        calc_timestep.cpp
        sim_control.cpp
        sim_control_MPI.cpp
        sim_init.cpp
        time_integrator.cpp
    )
    add_library(
        nested_sim_control

        sim_control_NG.cpp
        sim_control_NG_MPI.cpp
    )

else()
    add_library(
        sim_control

        calc_timestep.cpp
        sim_control.cpp
        sim_init.cpp
        time_integrator.cpp
    )
    add_library(
        nested_sim_control

        sim_control_NG.cpp
    )
endif(PION_PARALLEL)

target_link_libraries(nested_sim_control PUBLIC sim_control)
