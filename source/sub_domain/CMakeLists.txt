add_library(
    sub_domain 

    sub_domain_setup.cpp
    sub_domain_access.cpp
    sub_domain_comms.cpp
)

target_link_libraries(sub_domain PUBLIC ${MPI_CXX_LIBRARIES})
