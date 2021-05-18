
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(PION_TOOLS AND PION_PARALLEL)
  add_library(
    xray_emission

    xray_emission.cpp
  )

endif(PION_TOOLS AND PION_PARALLEL)
