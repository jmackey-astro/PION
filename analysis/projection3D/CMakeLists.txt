
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO AND PION_USE_FITS)
  add_executable(
    projection3D

    main_projection.cpp
    point_quantities.cpp
    sim_projection.cpp
  )

  target_compile_features(projection3D PRIVATE cxx_std_14)

  target_link_libraries(projection3D PRIVATE
    global
    dataio
    xray_emission
    image_io
  )
  if(PION_NESTED_GRID)
      target_link_libraries(projection3D PRIVATE nested_grid)
      target_compile_definitions(silocompare PRIVATE PION_NESTED)
  else()
      target_link_libraries(projection3D PRIVATE grid)
  endif(PION_NESTED_GRID)

  install( 
    TARGETS projection3D
    EXPORT projection3D-targets
    RUNTIME DESTINATION ${PION_INSTALL_BINDIR}
    LIBRARY DESTINATION ${PION_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${PION_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${PION_INSTALL_INCLUDEDIR}
  )

endif(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO AND PION_USE_FITS)
