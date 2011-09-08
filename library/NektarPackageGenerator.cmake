# CPack setup
IF (NOT ${NEKTAR_PACKAGE_GENERATOR} MATCHES "None")

	# Check if we have a valid Package generator selected
	IF (NOT ${NEKTAR_PACKAGE_GENERATOR} MATCHES "RPM"
			AND NOT ${NEKTAR_PACKAGE_GENERATOR} MATCHES "DEB")
		MESSAGE(SEND_ERROR "Unknown packaging generator: 
							${NEKTAR_PACKAGE_GENERATOR}")
	ENDIF()

	SET(CPACK_GENERATOR "${NEKTAR_PACKAGE_GENERATOR}")
	SET(CPACK_TOPLEVEL_TAG "${CMAKE_SYSTEM_PROCESSOR}")
	SET(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_PROCESSOR}")
	SET(CPACK_PACKAGE_NAME "nektar++")
	SET(CPACK_PACKAGE_VENDOR "Imperial College London")
	SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY 
	    "Spectral/hp Element Framework for solving PDEs")
	SET(CPACK_PACKAGE_INSTALL_DIRECTORY "nektar++")
	SET(CPACK_PACKAGE_CONTACT "Chris Cantwell <c.cantwell@imperial.ac.uk>")
	SET(CPACK_PACKAGE_VERSION ${NEKTAR_VERSION})
	SET(CPACK_PACKAGE_VERSION_MAJOR ${NEKTAR_VERSION_MAJOR})
	SET(CPACK_PACKAGE_VERSION_MINOR ${NEKTAR_VERSION_MINOR})
	SET(CPACK_PACKAGE_VERSION_PATCH ${NEKTAR_VERSION_PATCH})

	# RPM packaging
	IF (${NEKTAR_PACKAGE_GENERATOR} MATCHES "RPM")
		MESSAGE(STATUS "Generating Packaging for RPM")
		
		SET(CPACK_RPM_PACKAGE_URL "www.nektar.info")
		SET(CPACK_RPM_COMPONENT_INSTALL ON)
		SET(CPACK_RPM_PACKAGE_REQUIRES "fftw3, libboost_date_time1_44_0, libboost_filesystem1_44_0, libboost_iostreams1_44_0, libboost_system1_44_0, libboost_thread1_44_0, zlib")
		SET(CPACK_RPM_PACKAGE_DESCRIPTION "
The nektar++ packages provide a spectral/hp element framework for the numerical
solution of partial differential equations (PDEs). Demonstration codes are 
provided in the nektar++-demos package, along with a number of time-dependent 
solvers in the nektar++-solvers package.")
	ENDIF()

	# Finally, include the CPack module
	INCLUDE(CPack)

ENDIF()
