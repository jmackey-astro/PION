include(FindPackageHandleStandardArgs)

if(NOT SILO_FOUND)
    find_path( SILO_INCLUDE_DIR
        HINTS ./
            /usr/
            ${CMAKE_INSALL_PREFIX}/
            ${CMAKE_INSALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${SILO_DIR}/
        PATH_SUFFIXES include include/silo silo/include
        NAMES silo.h
        DOC "Silo headers"
    )
    find_library( SILO_LIBRARY
        HINTS ./
            /usr/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${SILO_DIR}/
        PATH_SUFFIXES lib lib64 silo/lib silo/lib64 lib/silo lib64/silo
        NAMES libsilo.a libsilo.la libsiloh5.so
        DOC "Silo libraries"
    )

    find_package_handle_standard_args(silo
            REQUIRED_VARS SILO_INCLUDE_DIR SILO_LIBRARY
    )            
endif(NOT SILO_FOUND)
