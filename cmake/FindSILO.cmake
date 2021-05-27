include(FindPackageHandleStandardArgs)

if(NOT SILO_FOUND)
    find_path( SILO_INCLUDE_DIR
        HINTS /
            /usr/
            ${CMAKE_INSALL_PREFIX}/
            ${CMAKE_INSALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${SILO_DIR}/
        PATH_SUFFIXES include
                    include/silo
                    include/x86_64-linux-gnu
                    include/x86_64-linux-gnu/silo
        NAMES silo.h
        DOC "Silo headers"
    )
    find_library( SILO_LIBRARY
        HINTS /
            /usr/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${SILO_DIR}/
        PATH_SUFFIXES lib lib64
                    lib/silo lib64/silo
                    lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
                    lib/x86_64-linux-gnu/silo lib64/x86_64-linux-gnu/silo
        NAMES libsilo.so libsilo.la libsilo.a libsiloh5.so
        DOC "Silo libraries"
    )

    find_package_handle_standard_args(SILO
            REQUIRED_VARS SILO_INCLUDE_DIR SILO_LIBRARY
    )            
endif(NOT SILO_FOUND)
