include(FindPackageHandleStandardArgs)

function(find_sundials_version)
    include(CheckCXXSymbolExists)
    include(CMakePushCheckState)
    set(CMAKE_REQUIRED_QUIET ON)

    CMAKE_PUSH_CHECK_STATE()

    set(CMAKE_REQUIRED_INCLUDES "${SUNDIALS_INCLUDE_DIR}")
    file(READ "${SUNDIALS_INCLUDE_DIR}/sundials_config.h" version_header)
    string(REGEX MATCH "#define SUNDIALS_VERSION \"[0-9]+.[0-9]+.[0-9]+\"" version_string_macrodef "${version_header}")
    string(REGEX MATCH "[0-9]+.[0-9]+.[0-9]+" SUNDIALS_VERSION_STRING "${version_string_macrodef}")
    string(REGEX MATCH "#define SUNDIALS_VERSION_MAJOR [0-9]+" version_major_macrodef "${version_header}")
    string(REGEX MATCH "[0-9]+" SUNDIALS_VERSION_MAJOR "${version_major_macrodef}")

    set(SUNDIALS_VERSION_MAJOR ${SUNDIALS_VERSION_MAJOR} PARENT_SCOPE)
    set(SUNDIALS_VERSION_STRING ${SUNDIALS_VERSION_STRING} PARENT_SCOPE)

    CMAKE_POP_CHECK_STATE()
endfunction()

if(NOT SUNDIALS_FOUND)
    find_path( SUNDIALS_INCLUDE_DIR
        HINTS /
            /usr/
            ${CMAKE_INSALL_PREFIX}/
            ${CMAKE_INSALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${SUNDIALS_DIR}/
        PATH_SUFFIXES include
                    include/sundials
                    include/x86_64-linux-gnu
                    include/x86_64-linux-gnu/sundials
        NAMES sundials_config.h
        DOC "Sundials headers"
    )
    find_sundials_version()    
    find_library( SUNDIALS_LIBRARY
        HINTS /
            /usr/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${SUNDIALS_DIR}/
        PATH_SUFFIXES lib lib64
                    lib/sundials lib64/sundials
                    lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
                    lib/x86_64-linux-gnu/sundials lib64/x86_64-linux-gnu/sundials
        NAMES libsundials_cvode.a libsundials_cvode.so
        DOC "Sundials libraries"
    )

    find_package_handle_standard_args(SUNDIALS
            FOUND_VAR SUNDIALS_FOUND
            REQUIRED_VARS SUNDIALS_INCLUDE_DIR SUNDIALS_LIBRARY SUNDIALS_VERSION_MAJOR
            VERSION_VAR SUNDIALS_VERSION_STRING
    )            
endif(NOT SUNDIALS_FOUND)
