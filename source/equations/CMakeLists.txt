add_library(
    equations

    eqns_base.cpp
    eqns_hydro_adiabatic.cpp
    eqns_mhd_adiabatic.cpp
)

target_link_libraries(
    equations
    PUBLIC 
    microphysics 
)

