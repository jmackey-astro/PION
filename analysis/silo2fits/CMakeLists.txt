
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(PION_TOOLS AND PION_PARALLEL AND PION_USE_FITS AND PION_USE_SILO)
  add_executable(
    silo2fits

    silo2fits.cpp
  )

  target_compile_features(silo2fits PUBLIC cxx_std_14)

  target_link_libraries(silo2fits PUBLIC 
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
  )

  install( 
    TARGETS silo2fits
    EXPORT silo2fits-targets
    RUNTIME DESTINATION ${PION_INSTALL_BINDIR}
    LIBRARY DESTINATION ${PION_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${PION_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${PION_INSTALL_INCLUDEDIR}
  )

endif(PION_TOOLS AND PION_PARALLEL AND PION_USE_FITS AND PION_USE_SILO)
