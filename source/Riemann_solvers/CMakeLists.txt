add_library(
    riemann_solvers

	findroot.cpp
	HLLD_MHD.cpp
	HLL_hydro.cpp
	riemann.cpp
	Riemann_FVS_hydro.cpp
	riemannMHD.cpp
	Roe_Hydro_ConservedVar_solver.cpp
	Roe_Hydro_PrimitiveVar_solver.cpp
	Roe_MHD_ConservedVar_solver.cpp
)

target_link_libraries(riemann_solvers PUBLIC equations)
target_link_libraries(riemann_solvers PRIVATE spdlog::spdlog)

