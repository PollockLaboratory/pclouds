cmake_minimum_required (VERSION 2.8.11 FATAL_ERROR)
project (pclouds CXX)
enable_testing (true)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

set (PCLOUDS_VERSION_MAJOR 0)
set (PCLOUDS_VERSION_MINOR 0)
set (PCLOUDS_VERSION_PATCH 0)
set (PCLOUDS_VERSION_TWEAK 20150601)
set (PCLOUDS_VERSION ${PCLOUDS_VERSION_MAJOR}.${PCLOUDS_VERSION_MINOR}.${PCLOUDS_VERSION_PATCH})
set (PCLOUDS_SOVERSION 0)

set (PCLOUDS_VENDOR "Pollock Laboratory")
set (PCLOUDS_DESCRIPTION "Identifying high copy repeat content in the genome using short k-mers.")
set (PCLOUDS_WEBSITE "http://pollocklab.github.io/pcloud")
set (PCLOUDS_SOURCE "https://github.com/pollocklab/pcloud")
set (PCLOUDS_AUTHOR_NAME "Daniel L. Polanco")
set (PCLOUDS_AUTHOR_EMAIL "daniel@blazingaddles.com")

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
