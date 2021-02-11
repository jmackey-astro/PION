add_library(
    comms

    comm_files.cpp
    comm_mpi.cpp
)

target_link_libraries(comms PUBLIC ${MPI_CXX_LIBRARIES})
