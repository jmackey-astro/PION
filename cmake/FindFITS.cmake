include(FindPackageHandleStandardArgs)

if(NOT FITS_FOUND)
    find_path( FITS_INCLUDE_DIR
        HINTS ./
            /usr/
            ${CMAKE_INSALL_PREFIX}/
            ${CMAKE_INSALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${FITS_DIR}/
        PATH_SUFFIXES include include/fits fits/include
        NAMES fitsio.h
        DOC "Fits headers"
    )
    find_library( FITS_LIBRARY
        HINTS ./
            /usr/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${FITS_DIR}/
        PATH_SUFFIXES lib lib64 fits/lib fits/lib64 lib/fits lib64/fits
        NAMES libcfitsio.a
        DOC "Fits libraries"
    )

    find_package_handle_standard_args(fits
            REQUIRED_VARS FITS_INCLUDE_DIR FITS_LIBRARY
    )            
endif(NOT FITS_FOUND)
