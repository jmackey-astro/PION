add_library(
    dataio

    dataio_base.cpp
    dataio_text.cpp
    file_status.cpp
    parameter_defs.cpp
    readparams.cpp
)

if(PION_USE_SILO)
    if(PION_PARALLEL)
        add_library(
            dataio_silo

            dataio_silo.cpp
            dataio_silo_MPI.cpp
            dataio_silo_utility.cpp
        )
    else()
        add_library(
            dataio_silo

            dataio_silo.cpp
            dataio_silo_utility.cpp
        )
    endif(PION_PARALLEL)

    target_link_libraries(dataio_silo PUBLIC ${SILO_LIBRARY})
    target_link_libraries(dataio PUBLIC dataio_silo)
endif(PION_USE_SILO)

if(PION_USE_FITS)
    if(PION_PARALLEL)
        add_library(
            dataio_fits

            dataio_fits.cpp
            dataio_fits_MPI.cpp
            utility_fits_class.cpp
        )
    else()
        add_library(
            dataio_fits

            dataio_fits.cpp
            utility_fits_class.cpp
        )
    endif(PION_PARALLEL)

    target_link_libraries(dataio_fits PUBLIC
        ${FITS_LIBRARY}
        ${ZLIB_LIBRARY}
        ${CURL_LIBRARY}
        ${BZIP2_LIBRARY}
        Threads::Threads
    )
    target_link_libraries(dataio PUBLIC dataio_fits)
endif(PION_USE_FITS)

target_link_libraries(dataio PUBLIC tools) 
