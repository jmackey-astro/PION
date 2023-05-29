include(FindPackageHandleStandardArgs)

if(NOT READLINE_FOUND)
    find_path( READLINE_INCLUDE_DIR
        HINTS /
            /usr/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${READLINE_DIR}/
            /usr/local/opt/readline/
        PATH_SUFFIXES include
                    include/readline
                    include/x86_64-linux-gnu
                    include/x86_64-linux-gnu/readline
        NAMES readline.h
        DOC "Readline headers"
    )
    find_library( READLINE_LIBRARY
        HINTS /
            /usr/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${READLINE_DIR}/
            /usr/local/opt/readline/
        PATH_SUFFIXES lib lib64
                    lib/readline lib64/readline
                    lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
                    lib/x86_64-linux-gnu/readline lib64/x86_64-linux-gnu/readline
        NAMES libreadline.so libreadline.la libreadline.a
        DOC "Readline libraries"
    )

    find_package_handle_standard_args(READLINE
            REQUIRED_VARS READLINE_INCLUDE_DIR READLINE_LIBRARY
    )            
endif(NOT READLINE_FOUND)
