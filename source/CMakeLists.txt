# add internal libraries
add_subdirectory(boundaries)
add_subdirectory(comms)
#add_subdirectory(contrib)
add_subdirectory(coord_sys)
add_subdirectory(dataIO)
add_subdirectory(decomposition)
add_subdirectory(equations)
add_subdirectory(ics)
add_subdirectory(microphysics)
add_subdirectory(raytracing)
add_subdirectory(grid)
add_subdirectory(Riemann_solvers)
add_subdirectory(sim_control)
add_subdirectory(spatial_solvers)
add_subdirectory(tools)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


if(PION_NESTED_GRID)
    add_executable(
        pion-ng

        main_NG_MPI.cpp
        sim_params.cpp
        constants.cpp
    ) 

    add_executable(
        icgen-ng

        main_icgen_NG_MPI.cpp
        sim_params.cpp
        constants.cpp
    ) 

    target_compile_features(pion-ng PUBLIC cxx_std_14)
    target_compile_features(icgen-ng PUBLIC cxx_std_14)

    # pion-ng internal dependencies
    target_link_libraries(pion-ng PUBLIC 
        sim_control
        decomposition
        raytracing
    )

    # icgen-ng internal dependencies
    target_link_libraries(icgen-ng PUBLIC 
        ics
        decomposition
        raytracing
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
    add_executable(
        pion-ug

        mainMPI.cpp
        sim_params.cpp
        constants.cpp
    ) 

    add_executable(
        icgen-ug

        main_icgen_parallel.cpp
        sim_params.cpp
        constants.cpp
    ) 

    target_compile_features(pion-ug PUBLIC cxx_std_14)
    target_compile_features(icgen-ug PUBLIC cxx_std_14)

    # pion-ug internal dependencies
    target_link_libraries(pion-ug PUBLIC 
        sim_control
        decomposition
        raytracing
    )

    # icgen-ug internal dependencies
    target_link_libraries(icgen-ug PUBLIC 
        ics
        decomposition
        raytracing
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
