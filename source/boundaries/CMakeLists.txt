if(PION_PARALLEL)
    add_library(
        boundaries

        assign_update_bcs.cpp
        assign_update_bcs_MPI.cpp
        axisymmetric_boundaries.cpp
        double_Mach_ref_boundaries.cpp
        fixed_boundaries.cpp
        inflow_boundaries.cpp
        jet_boundaries.cpp
        jetreflect_boundaries.cpp
        MCMD_boundaries.cpp
        oneway_out_boundaries.cpp
        outflow_boundaries.cpp
        periodic_boundaries.cpp
        periodic_boundaries_MPI.cpp
        reflecting_boundaries.cpp
        RT_MPI_boundaries.cpp
        stellar_wind_boundaries.cpp
    )
    add_library(
        nested_boundaries

        assign_update_bcs_NG.cpp
        assign_update_bcs_NG_MPI.cpp
        NG_BC89flux.cpp
        NG_coarse_to_fine_boundaries.cpp
        NG_fine_to_coarse_boundaries.cpp
        NG_MPI_BC89flux.cpp
        NG_MPI_coarse_to_fine_boundaries.cpp
        NG_MPI_fine_to_coarse_boundaries.cpp
    )
else()
    add_library(
        boundaries

        assign_update_bcs.cpp
        axisymmetric_boundaries.cpp
        double_Mach_ref_boundaries.cpp
        fixed_boundaries.cpp
        inflow_boundaries.cpp
        jet_boundaries.cpp
        jetreflect_boundaries.cpp
        oneway_out_boundaries.cpp
        outflow_boundaries.cpp
        periodic_boundaries.cpp
        reflecting_boundaries.cpp
        stellar_wind_boundaries.cpp
    )
    add_library(
        nested_boundaries

        assign_update_bcs_NG.cpp
        NG_BC89flux.cpp
        NG_coarse_to_fine_boundaries.cpp
        NG_fine_to_coarse_boundaries.cpp
    )
endif(PION_PARALLEL)

target_link_libraries(nested_boundaries PUBLIC
    boundaries
    nested_grid
)

target_link_libraries(boundaries PUBLIC grid)
