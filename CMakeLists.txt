cmake_minimum_required (VERSION 3.0 FATAL_ERROR)


## project overview
project (pclouds CXX C)
set (PCLOUDS_VERSION "0.0.0")
set (PCLOUDS_DESCRIPTION "Identification and annotation of high copy repeat content in the genome.")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

## options
set (PCLOUDS_TESTING FALSE "Compile tests for sequence.")

## prerequisites

## build
include_directories (include)
add_library (sequence
	sources/annotate.cpp)

if (PCLOUDS_TEST)
	enable_testing(true)
endif (PCLOUDS_TEST)
