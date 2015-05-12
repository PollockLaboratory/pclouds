include (CMakePackageConfigHelpers)
write_basic_package_version_file (
	${CMAKE_CURRENT_BINARY_DIR}/cmake/PCLOUDSConfigVersion.cmake
	VERSION ${PCLOUDS_VERSION}
	COMPATIBILITY AnyNewerVersion)

export (EXPORT PCLOUDSTargets
	FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/PCLOUDSTargets.cmake)

configure_file (
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/PCLOUDSConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cmake/PCLOUDSConfig.cmake
	COPYONLY)

install (FILES
	${CMAKE_CURRENT_BINARY_DIR}/cmake/PCLOUDSConfigVersion.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cmake/PCLOUDSTargets.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cmake/PCLOUDSConfig.cmake
	DESTINATION ${INSTALL_PACKAGE_DIR}
	COMPONENT devel)
