add_library(
    tools

    command_line_interface.cpp
    interpolate.cpp
    mem_manage.cpp
    reporting.cpp
    timer.cpp
)

target_link_libraries(tools PUBLIC GSL::gsl)

if(PION_PARALLEL)
    target_link_libraries(tools PUBLIC comms)
endif(PION_PARALLEL)
