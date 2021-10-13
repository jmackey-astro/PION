
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO)
  add_executable(
    silo2text

    silo2text.cpp
  )

  target_compile_features(silo2text PUBLIC cxx_std_14)

  target_link_libraries(silo2text PUBLIC 
    global
    dataio
    nested_grid
  )

  install( 
    TARGETS silo2text
    EXPORT silo2text-targets
    RUNTIME DESTINATION ${PION_INSTALL_BINDIR}
    LIBRARY DESTINATION ${PION_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${PION_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${PION_INSTALL_INCLUDEDIR}
  )

endif(PION_TOOLS AND PION_PARALLEL AND PION_USE_SILO)
