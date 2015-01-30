cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

add_executable (pclouds PClouds.cpp
	include/macrodefine.h
	#cloudbuilder/build_clouds.cpp
	cloudbuilder/build_clouds_recursively.cpp
	filereader/readfile.h
	filereader/readfile.cpp
	genomeannotator/annotate_genome.cpp
	stringhandler/stringhandle.h
	stringhandler/stringhandle.cpp
	)
