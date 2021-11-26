# add internal libraries
add_subdirectory(boundaries)
#add_subdirectory(contrib)
add_subdirectory(coord_sys)
add_subdirectory(dataIO)
add_subdirectory(equations)
add_subdirectory(grid)
add_subdirectory(ics)
add_subdirectory(microphysics)
add_subdirectory(raytracing)
add_subdirectory(Riemann_solvers)
add_subdirectory(sim_control)
add_subdirectory(spatial_solvers)
add_subdirectory(tools)

if (PION_PARALLEL)
    add_subdirectory(sub_domain)
endif(PION_PARALLEL)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(
  global

  sim_params.cpp
  constants.cpp
) 

target_link_libraries(global INTERFACE dataio tools)# spdlog::spdlog)

if(PION_OMP)
    target_link_libraries(global INTERFACE OpenMP::OpenMP_CXX)
endif(PION_OMP)

if(PION_NESTED_GRID)
    add_executable(pion-ng main.cpp) 
    add_executable(icgen-ng main_icgen.cpp)

    target_link_libraries(pion-ng PRIVATE global)
    target_link_libraries(icgen-ng PRIVATE global)

    target_compile_features(pion-ng PRIVATE cxx_std_14)
    target_compile_features(icgen-ng PRIVATE cxx_std_14)

    target_compile_definitions(pion-ng PRIVATE PION_NESTED)
    target_compile_definitions(icgen-ng PRIVATE PION_NESTED)

    # pion-ng internal dependencies
    target_link_libraries(pion-ng PRIVATE
        nested_sim_control
        nested_grid
    )

    # icgen-ng internal dependencies
    target_link_libraries(icgen-ng PRIVATE
        ics
        nested_grid
    )

    #
    # Set up installation directories
    #
    install( 
        TARGETS pion-ng icgen-ng
        EXPORT pion-ng-targets
        RUNTIME DESTINATION ${PION_INSTALL_BINDIR}
        LIBRARY DESTINATION ${PION_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${PION_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${PION_INSTALL_INCLUDEDIR}
    )
endif(PION_NESTED_GRID)

if(PION_UNIFORM_GRID)
    add_executable(pion-ug main.cpp)
    add_executable(icgen-ug main_icgen.cpp) 

    target_link_libraries(pion-ug PRIVATE global)
    target_link_libraries(icgen-ug PRIVATE global)

    target_compile_features(pion-ug PRIVATE cxx_std_14)
    target_compile_features(icgen-ug PRIVATE cxx_std_14)

    # pion-ug internal dependencies
    target_link_libraries(pion-ug PRIVATE 
        sim_control
        grid 
    )

    # icgen-ug internal dependencies
    target_link_libraries(icgen-ug PRIVATE 
        ics
        grid
    )

    #
    # Set up installation directories
    #
    install( 
        TARGETS pion-ug icgen-ug
        EXPORT pion-ug-targets
        RUNTIME DESTINATION ${PION_INSTALL_BINDIR}
        LIBRARY DESTINATION ${PION_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${PION_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${PION_INSTALL_INCLUDEDIR}
    )
endif(PION_UNIFORM_GRID)
