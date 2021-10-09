# add internal libraries
add_subdirectory(boundaries)
add_subdirectory(comms)
#add_subdirectory(contrib)
add_subdirectory(coord_sys)
add_subdirectory(dataIO)
add_subdirectory(equations)
add_subdirectory(ics)
add_subdirectory(microphysics)
add_subdirectory(raytracing)
add_subdirectory(grid)
add_subdirectory(Riemann_solvers)
add_subdirectory(sim_control)
add_subdirectory(spatial_solvers)
add_subdirectory(tools)
if (PION_PARALLEL)
    add_subdirectory(decomposition)
endif(PION_PARALLEL)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


if(PION_NESTED_GRID)
    if(PION_PARALLEL)
        add_executable(
            pion-ng

            main_NG.cpp
            sim_params.cpp
            constants.cpp
        ) 
        add_executable(
            icgen-ng

            main_icgen_NG.cpp
            sim_params.cpp
            constants.cpp
        ) 
        target_link_libraries(pion-ng PUBLIC decomposition)
        target_link_libraries(icgen-ng PUBLIC decomposition)
    else()
        add_executable(
            pion-ng

            main_NG.cpp
            sim_params.cpp
            constants.cpp
        ) 
        add_executable(
            icgen-ng

            main_icgen_NG.cpp
            sim_params.cpp
            constants.cpp
        )
    endif(PION_PARALLEL)

    target_compile_features(pion-ng PUBLIC cxx_std_14)
    target_compile_features(icgen-ng PUBLIC cxx_std_14)

    # pion-ng internal dependencies
    target_link_libraries(pion-ng PUBLIC 
        nested_sim_control
        raytracing
        nested_boundaries
    )

    # icgen-ng internal dependencies
    target_link_libraries(icgen-ng PUBLIC 
        ics
        raytracing
        nested_boundaries
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
    if(PION_PARALLEL)
        add_executable(
            pion-ug

            main.cpp
            sim_params.cpp
            constants.cpp
        )
        add_executable(
            icgen-ug

            main_icgen.cpp
            sim_params.cpp
            constants.cpp
        ) 
        target_link_libraries(pion-ug PUBLIC decomposition)
        target_link_libraries(icgen-ug PUBLIC decomposition)
    else()
        add_executable(
            pion-ug

            main.cpp
            sim_params.cpp
            constants.cpp
        )
        add_executable(
            icgen-ug

            main_icgen.cpp
            sim_params.cpp
            constants.cpp
        ) 
    endif(PION_PARALLEL)

    target_compile_features(pion-ug PUBLIC cxx_std_14)
    target_compile_features(icgen-ug PUBLIC cxx_std_14)

    # pion-ug internal dependencies
    target_link_libraries(pion-ug PUBLIC 
        sim_control
        raytracing
        boundaries
    )

    # icgen-ug internal dependencies
    target_link_libraries(icgen-ug PUBLIC 
        ics
        raytracing
        boundaries
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

add_library(
  global

  sim_params.cpp
  constants.cpp
)

#if(PION_OMP)
#    target_link_libraries(global PUBLIC OpenMP::OpenMP_CXX)
#endif(PION_OMP)

