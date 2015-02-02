cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -Wall")

add_executable (pclouds pclouds.cpp
	include/macrodefine.h
	cloudbuilder/build_clouds.cpp
	filereader/readfile.h
	filereader/readfile.cpp
	genomeannotator/annotate_genome.cpp
	stringhandler/stringhandle.h
	stringhandler/stringhandle.cpp
	)

add_subdirectory (processors)
