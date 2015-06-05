# Packaging
SET(NEKTAR_PACKAGE_GENERATOR "None" CACHE STRING 
    "Support Packaging: RPM or DEB" )
MARK_AS_ADVANCED(NEKTAR_PACKAGE_GENERATOR)

# Define components
SET(CPACK_COMPONENTS_ALL dev lib demos extra-demos)

#SET(CPACK_RPM_PACKAGE_DEBUG ON)
SET(CPACK_COMPONENT_LIB_NAME "nektar++-lib")
SET(CPACK_COMPONENT_LIB_DISPLAY_NAME "nektar++-lib")
SET(CPACK_COMPONENT_LIB_DESCRIPTION "Nektar++ Libraries")

SET(CPACK_COMPONENT_DEV_NAME "nektar++-dev")
SET(CPACK_COMPONENT_DEV_DISPLAY_NAME "nektar++-dev")
SET(CPACK_COMPONENT_DEV_DESCRIPTION "Development files for Nektar++")
SET(CPACK_COMPONENT_DEV_DEPENDS lib)

SET(CPACK_COMPONENT_DEMOS_NAME "nektar++-demos")
SET(CPACK_COMPONENT_DEMOS_DISPLAY_NAME "nektar++-demos")
SET(CPACK_COMPONENT_DEMOS_DESCRIPTION "Framework demonstration binaries")
SET(CPACK_COMPONENT_DEMOS_DEPENDS lib)

SET(CPACK_COMPONENT_EXTRA-DEMOS_NAME "nektar++-extra-demos")
SET(CPACK_COMPONENT_EXTRA-DEMOS_DISPLAY_NAME "nektar++-extra-demos")
SET(CPACK_COMPONENT_EXTRA-DEMOS_DESCRIPTION 
	"Framework extra demonstration binaries")
SET(CPACK_COMPONENT_EXTRA-DEMOS_DEPENDS lib)

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
	SET(CPACK_PACKAGE_DESCRIPTION "
The nektar++ packages provide a spectral/hp element framework for the numerical
solution of partial differential equations (PDEs). Demonstration codes are 
provided in the nektar++-demos package, along with a number of time-dependent 
solvers in the nektar++-solvers package.")
	SET(CPACK_PACKAGE_INSTALL_DIRECTORY "nektar++")
	SET(CPACK_PACKAGE_CONTACT "Chris Cantwell <c.cantwell@imperial.ac.uk>")
	SET(CPACK_PACKAGE_VERSION ${NEKTAR_VERSION})
	SET(CPACK_PACKAGE_VERSION_MAJOR ${NEKTAR_VERSION_MAJOR})
	SET(CPACK_PACKAGE_VERSION_MINOR ${NEKTAR_VERSION_MINOR})
	SET(CPACK_PACKAGE_VERSION_PATCH ${NEKTAR_VERSION_PATCH})
    SET(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/LICENSE)

	# RPM packaging
	IF (${NEKTAR_PACKAGE_GENERATOR} MATCHES "RPM")
		MESSAGE(STATUS "Generating Packaging for RPM")
		
		SET(CPACK_RPM_PACKAGE_URL "www.nektar.info")
		SET(CPACK_RPM_COMPONENT_INSTALL ON)
		SET(CPACK_RPM_PACKAGE_REQUIRES "fftw3, libboost_date_time1_44_0, libboost_filesystem1_44_0, libboost_iostreams1_44_0, libboost_system1_44_0, libboost_thread1_44_0, libboost_timer1_44_0, zlib")
		SET(CPACK_RPM_PACKAGE_DESCRIPTION "
The nektar++ packages provide a spectral/hp element framework for the numerical
solution of partial differential equations (PDEs). Demonstration codes are 
provided in the nektar++-demos package, along with a number of time-dependent 
solvers in the nektar++-solvers package.")
	ENDIF()

    IF (${NEKTAR_PACKAGE_GENERATOR} MATCHES "DEB")
        MESSAGE(STATUS "Generating Packaging for DEB")
        SET(CPACK_DEB_PACKAGE_URL "www.nektar.info")
        SET(CPACK_DEB_COMPONENT_INSTALL ON)
        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libfftw3-3,libboost-date-time1.42.0,libboost-filesystem1.42.0,libboost-iostreams1.42.0,libboost-program-options1.42.0,libboost-system1.42.0,libboost-thread1.42.0,libboost-timer1.42.0,zlib1g")
        SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION 
            "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}
        ${CPACK_PACKAGE_DESCRIPTION}")
        SET(CPACK_DEBIAN_PACKAGE_SECTION "devel")
        SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
	ENDIF()
        
    # Finally, include the CPack module
	INCLUDE(CPack)

ENDIF()
