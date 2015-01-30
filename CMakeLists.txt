cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

add_subdirectory (cloudbuilder)
add_subdirectory (filereader)
add_subdirectory (genomeannotator)
add_subdirectory (processors)
add_subdirectory (stringhandler)

add_executable (pclouds PClouds.cpp)
