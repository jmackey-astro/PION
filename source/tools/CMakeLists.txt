add_library( tools
  
  command_line_interface.cpp
  interpolate.cpp
  mem_manage.cpp
  reporting.cpp
  timer.cpp
)

#target_link_libraries(tools PUBLIC GSL::gsl)
target_link_libraries(tools PRIVATE Boost::boost)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Curses_FOUND: ${Curses_FOUND}")
    if (Curses_FOUND)
        target_link_libraries(tools PRIVATE ${CURSES_LIBRARIES})
    endif(Curses_FOUND)
    target_link_libraries(tools PRIVATE ${READLINE_LIBRARY})
endif(CMAKE_BUILD_TYPE MATCHES "Debug")
