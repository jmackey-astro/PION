macro(custom_add_test)
    set(OPTIONS EXCLUDE)
    set(SINGLE_VALUED NAME SCRIPT ARGUMENTS)
    cmake_parse_arguments(TEST "${OPTIONS}" "${SINGLE_VALUED}" "" ${ARGN})

    add_custom_target(
        test-${TEST_NAME}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_SCRIPT} ${TEST_ARGUMENTS}
        VERBATIM
    )
    if (NOT ${TEST_EXCLUDE})
        add_test(
            NAME ${TEST_NAME}
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_SCRIPT} ${TEST_ARGUMENTS}
        )
    endif (NOT ${TEST_EXCLUDE})
endmacro(custom_add_test)

if (PION_PARALLEL)
    if (PION_UNIFORM_GRID)
        custom_add_test(
            NAME uniform
            SCRIPT dmr_ug.sh
        )
        custom_add_test(
            NAME pbc
            SCRIPT otv_pbc.sh
        )
        custom_add_test(
            NAME 3d-uniform
            SCRIPT blastwave_3d_ug.sh
        )
        custom_add_test(
            NAME 3d-uniform-pbc
            SCRIPT blastwave_3d_ug_pbc.sh
        )
        custom_add_test(
            NAME dte-mpi
            SCRIPT run-2d-rhd.sh
            EXCLUDE
        )
        custom_add_test(
            NAME dte-omp
            SCRIPT run-2d-rhd.sh
            ARGUMENTS 1
            EXCLUDE
        )
        custom_add_test(
            NAME dte-hybrid
            SCRIPT run-2d-rhd.sh
            ARGUMENTS 2
        )
        custom_add_test(
            NAME binary-3d-ug
            SCRIPT binary_3d_ug.sh
            EXCLUDE
        )
        custom_add_test(
            NAME binary-3d-ug-adi
            SCRIPT binary_3d_adi_ug.sh
            EXCLUDE
        )
    endif(PION_UNIFORM_GRID)
    if (PION_NESTED_GRID)
        custom_add_test(
            NAME nested
            SCRIPT mhd_blastwave_ng.sh
        )
        custom_add_test(
            NAME 3d-nested
            SCRIPT blastwave_3d_ng.sh
        )
        custom_add_test(
            NAME 3d-nested-pbc
            SCRIPT blastwave_3d_ng_pbc.sh
        )
        custom_add_test(
            NAME binary-3d-ng
            SCRIPT binary_3d_ng.sh
            EXCLUDE
        )
        custom_add_test(
            NAME binary-3d-ng-adi
            SCRIPT binary_3d_adi_ng.sh
            EXCLUDE
        )
    endif(PION_NESTED_GRID)
endif(PION_PARALLEL)
