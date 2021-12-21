set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO)
  add_executable(
    silocompare

    silocompare.cpp
  )

  target_compile_features(silocompare PRIVATE cxx_std_14)

  target_link_libraries(silocompare PRIVATE
    global
    dataio
    spdlog::spdlog
  )

  if (PION_NESTED_GRID)
      target_link_libraries(silocompare PRIVATE nested_grid)
      target_compile_definitions(silocompare PRIVATE PION_NESTED)
  else()
      target_link_libraries(silocompare PRIVATE grid)
  endif(PION_NESTED_GRID)

  install( 
    TARGETS silocompare
    EXPORT silocompare-targets
    RUNTIME DESTINATION ${PION_INSTALL_BINDIR}
    LIBRARY DESTINATION ${PION_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${PION_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${PION_INSTALL_INCLUDEDIR}
  )

endif(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO)
