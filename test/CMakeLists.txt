add_custom_target(ug-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/dmr_ug.sh
    VERBATIM
)

add_custom_target(ng-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/mhd_blastwave_ng.sh
    VERBATIM
)

add_custom_target(pbc-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/otv_pbc.sh
    VERBATIM
)

add_custom_target(3d-ug-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/blastwave_3d_ug.sh
    VERBATIM
)

add_custom_target(3d-ng-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/blastwave_3d_ng.sh
    VERBATIM
)
add_custom_target(3d-ug-pbc-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/blastwave_3d_ug_pbc.sh
    VERBATIM
)

add_custom_target(3d-ng-pbc-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/blastwave_3d_ng_pbc.sh
    VERBATIM
)


add_custom_target(dte_test_mpi
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/run-2d-rhd.sh
    VERBATIM
)

add_custom_target(dte_test_omp
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/run-2d-rhd.sh 1
    VERBATIM
)

add_custom_target(dte_test_hybrid
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/run-2d-rhd.sh 2
    VERBATIM
)

add_custom_target(binary-3d-ug-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/binary_3d_ug.sh
    VERBATIM
)

add_custom_target(binary-3d-ng-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/binary_3d_ng.sh
    VERBATIM
)
add_custom_target(binary-3d-ug-adi-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/binary_3d_adi_ug.sh
    VERBATIM
)

add_custom_target(binary-3d-ng-adi-test
    COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/binary_3d_adi_ng.sh
    VERBATIM
)




