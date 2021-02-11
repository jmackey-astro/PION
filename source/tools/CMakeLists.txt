add_library(
    tools

    command_line_interface.cpp
    interpolate.cpp
    mem_manage.cpp
    reporting.cpp
    timer.cpp
)

target_link_libraries(tools PUBLIC 
    GSL::gsl
    comms
)
