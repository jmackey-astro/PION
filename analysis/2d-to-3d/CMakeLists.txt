set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO AND PION_NESTED_GRID)
  add_executable(
    map-2d-to-3d

    map-2d-to-3d.cpp
  )

  target_compile_features(map-2d-to-3d PRIVATE cxx_std_14)

  target_link_libraries(map-2d-to-3d PRIVATE
    global
    dataio
    spdlog::spdlog
  )

  if (PION_NESTED_GRID)
      target_link_libraries(map-2d-to-3d PRIVATE nested_grid)
      target_compile_definitions(map-2d-to-3d PRIVATE PION_NESTED)
  else()
      target_link_libraries(map-2d-to-3d PRIVATE grid)
  endif(PION_NESTED_GRID)

  install( 
    TARGETS map-2d-to-3d
    EXPORT map-2d-to-3d-targets
    RUNTIME DESTINATION ${PION_INSTALL_BINDIR}
    LIBRARY DESTINATION ${PION_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${PION_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${PION_INSTALL_INCLUDEDIR}
  )

endif(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO AND PION_NESTED_GRID)
