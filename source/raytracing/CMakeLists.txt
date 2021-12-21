if(PION_PARALLEL)
    add_library(
        raytracing

        raytracer_SC.cpp
        raytracer_SC_pllel.cpp
    )
else()
    add_library(
        raytracing

        raytracer_SC.cpp
    )
endif(PION_PARALLEL)

target_link_libraries(raytracing PRIVATE global boundaries)
target_link_libraries(raytracing PRIVATE spdlog::spdlog)

