cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

set (PROJECT_NAME pclouds)
project (${PROJECT_NAME})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

option (test "Build all tests." OFF)

add_library (cloudbuilder_lib cloudbuilder/build_clouds.hpp cloudbuilder/build_clouds.cpp)
add_library (filereader_lib filereader/readfile.hpp filereader/readfile.cpp)
add_library (genomeannotator_lib genomeannotator/annotate_genome.cpp)
add_library (stringhandler_lib stringhandler/stringhandle.hpp stringhandler/stringhandle.cpp)

add_executable (pclouds macrodefine.h pclouds/pclouds.cpp)
target_link_libraries (pclouds cloudbuilder_lib)
target_link_libraries (pclouds filereader_lib)
target_link_libraries (pclouds genomeannotator_lib)
target_link_libraries (pclouds stringhandler_lib)

add_subdirectory (processors)

if (test)
	enable_testing ()

	add_executable (cloudbuilder cloudbuilder/test.cpp)
	target_link_libraries (cloudbuilder cloudbuilder_lib)
	target_link_libraries (cloudbuilder filereader_lib)
	target_link_libraries (cloudbuilder genomeannotator_lib)
	target_link_libraries (cloudbuilder stringhandler_lib)

	add_test (cloud_builder_test cloudbuilder)
endif ()
