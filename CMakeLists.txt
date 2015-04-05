cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

set (PROJECT_NAME pclouds)
project (${PROJECT_NAME})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

option (PCLOUDS_TEST "Build all tests." OFF)

add_library (cloudbuilder_lib cloudbuilder/build_clouds.hpp cloudbuilder/build_clouds.cpp)
add_library (filereader_lib filereader/readfile.hpp filereader/readfile.cpp)
add_library (genomeannotator_lib genomeannotator/annotate_genome.cpp)
add_library (stringhandler_lib stringhandler/stringhandle.hpp stringhandler/stringhandle.cpp)

add_executable (pclouds macrodefine.h pclouds/pclouds.cpp)
target_link_libraries (pclouds
	cloudbuilder_lib
	filereader_lib
	genomeannotator_lib
	stringhandler_lib)

add_subdirectory (processors)

if (PCLOUDS_TEST)
	enable_testing ()
	file (COPY controlfile DESTINATION ./)
	#add_executable (cloudbuilder cloudbuilder/test.cpp
	#	cloudbuilder_lib
	#	filereader_lib
	#	genomeannotator_lib
	#	stringhandler_lib)
	#add_test (cloud_builder_test cloudbuilder)
endif (PCLOUDS_TEST)
