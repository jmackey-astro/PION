add_library(
    ics

	basic_tests.cpp                          
	blast_wave.cpp                           
	icgen_base.cpp                           
	inside_sphere.cpp                        
	jet.cpp                                  
	laser_ablation.cpp                       
	photoevaporating_clump.cpp               
	photoevaporating_multiclumps.cpp         
	photoevaporating_random_clumps.cpp       
	radiative_shock.cpp                      
	read_1Dto2D.cpp                          
	shock_cloud.cpp                          
	shock_tube.cpp                           
	spherical_clump.cpp                      
	StarBench_test.cpp                       
)

target_link_libraries(ics PRIVATE spdlog::spdlog)
