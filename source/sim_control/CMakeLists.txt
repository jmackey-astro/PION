add_library(
    sim_control

    calc_timestep.cpp
    sim_control.cpp
    sim_control_MPI.cpp
    sim_control_NG.cpp
    sim_control_NG_MPI.cpp
    sim_init.cpp
    time_integrator.cpp
)
