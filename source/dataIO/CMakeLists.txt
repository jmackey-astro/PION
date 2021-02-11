add_library(
    dataio

    dataio_base.cpp
    dataio_fits.cpp
    dataio_fits_MPI.cpp
    dataio_silo.cpp
    dataio_silo_MPI.cpp
    dataio_silo_utility.cpp
    dataio_text.cpp
    file_status.cpp
    parameter_defs.cpp
    readparams.cpp
    utility_fits_class.cpp
)

target_link_libraries(dataio PUBLIC 
    ${SILO_LIBRARY}
    tools
) 
