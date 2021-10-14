add_library(
    microphysics

	cooling.cpp
	cooling_SD93_cie.cpp
	cvode_integrator.cpp
	hydrogen_mp.cpp
	hydrogen_photoion.cpp
	hydrogen_recomb_Hummer94.cpp
	integrator.cpp
	microphysics_base.cpp
	mp_only_cooling.cpp
	MPv0.cpp
	MPv1.cpp
	MPv2.cpp
	MPv3.cpp
	MPv4.cpp
	MPv5.cpp
	MPv6.cpp
	MPv7.cpp
	MPv8.cpp
	MPv9.cpp
)

target_link_libraries(microphysics PUBLIC raytracing)
target_link_libraries(microphysics PRIVATE tools ${SUNDIALS_LIBRARY})
