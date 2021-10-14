set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO AND PION_USE_FITS)
  add_executable(
    projection2D

    angle_projection.cpp
    image_io.cpp
    perp_projection.cpp
    project2D.cpp
  )

  target_compile_features(projection2D PRIVATE cxx_std_14)

  target_link_libraries(projection2D PRIVATE
    global
    dataio
    xray_emission
  )
  if(PION_NESTED_GRID)
      target_link_libraries(projection2D PRIVATE nested_grid)
      target_compile_definitions(silocompare PRIVATE PION_NESTED)
  else()
      target_link_libraries(projection2D PRIVATE grid)
  endif()

  install( 
    TARGETS projection2D
    EXPORT projection2D-targets
    RUNTIME DESTINATION ${PION_INSTALL_BINDIR}
    LIBRARY DESTINATION ${PION_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${PION_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${PION_INSTALL_INCLUDEDIR}
  )

  add_library(
    image_io

    image_io.cpp
  )
  target_link_libraries(image_io PRIVATE dataio)

endif(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO AND PION_USE_FITS)
