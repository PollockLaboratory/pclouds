cmake_minimum_required (VERSION 3.0 FATAL_ERROR)

set (PROJECT_NAME pclouds)
project (${PROJECT_NAME})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

option (PCLOUDS_TEST "Build all tests." OFF)

if (PCLOUDS_TEST)
endif (PCLOUDS_TEST)
