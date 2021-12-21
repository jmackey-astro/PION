add_library(
    coord_sys

    VectorOps.cpp
    VectorOps_spherical.cpp
)

target_link_libraries(coord_sys PRIVATE spdlog::spdlog)
