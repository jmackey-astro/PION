include(FindPackageHandleStandardArgs)

if(NOT FITS_FOUND)
    find_path( FITS_INCLUDE_DIR
        HINTS .
            /usr/
            ${CMAKE_INSALL_PREFIX}/
            ${CMAKE_INSALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${FITS_DIR}/
        PATH_SUFFIXES include
                    include/fits
                    include/x86_64-linux-gnu
                    include/x86_64-linux-gnu/fits
        NAMES fitsio.h
        DOC "Fits headers"
    )
    find_library( FITS_LIBRARY
        HINTS .
            /usr/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${FITS_DIR}/
        PATH_SUFFIXES lib lib64
                    lib/fits lib64/fits
                    lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
                    lib/x86_64-linux-gnu/fits lib64/x86_64-linux-gnu/fits
        NAMES libcfitsio.a
        DOC "Fits libraries"
    )

    find_package_handle_standard_args(fits
            REQUIRED_VARS FITS_INCLUDE_DIR FITS_LIBRARY
    )            
endif(NOT FITS_FOUND)
