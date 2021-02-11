add_library(
    raytracing

    raytracer_SC.cpp
    raytracer_SC_pllel.cpp
)

target_link_libraries(raytracing PUBLIC boundaries)
