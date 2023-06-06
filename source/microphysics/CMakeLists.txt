
if(PION_MPV10)
  add_library(
    mpv10

    # The following lines are added by Arun Mathew
    MPv10.cpp
    mellema_cooling.cpp
    atomic_physics_database.cpp
    coll_ionise_recomb.cpp
    chianti_cooling.cpp
  )
endif(PION_MPV10)

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

if(PION_MPV10)
  target_link_libraries(microphysics PUBLIC mpv10)
endif(PION_MPV10)
target_link_libraries(microphysics PUBLIC raytracing)
target_link_libraries(microphysics PRIVATE tools ${SUNDIALS_LIBRARY})
target_link_libraries(microphysics PRIVATE spdlog::spdlog)



