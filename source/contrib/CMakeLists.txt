add_library(
    contrib

    HD_2D_ShockCloud.cpp
    HD_3D_ShockCloud.cpp
    read_BBurkhart_data.cpp
)

target_link_libraries(contrib PRIVATE spdlog::spdlog)
