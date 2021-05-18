
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO AND PION_USE_FITS)
  add_executable(
    projection3D

    main_projection.cpp
    point_quantities.cpp
    sim_projection.cpp
  )

  target_compile_features(projection3D PUBLIC cxx_std_14)

  target_link_libraries(projection3D PUBLIC 
    global
    dataio
    dataio_silo
    dataio_fits
    grid
    nested_grid
    decomposition
    microphysics
    coord_sys
    equations
    spatial_solvers
    raytracing
    boundaries
    nested_boundaries
    xray_emission
    image_io
  )

  install( 
    TARGETS projection3D
    EXPORT projection3D-targets
    RUNTIME DESTINATION ${PION_INSTALL_BINDIR}
    LIBRARY DESTINATION ${PION_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${PION_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${PION_INSTALL_INCLUDEDIR}
  )

endif(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO AND PION_USE_FITS)
