########################################################################
#
# ThirdParty configuration for Nektar++
#
# OpenCascade
#
########################################################################

IF(NEKTAR_USE_MESHGEN)

    OPTION(NEKTAR_USE_CFI "Use CFI as cad engine." OFF)

    IF (NEKTAR_USE_CFI)

	SET(TEST_ENV $ENV{FEGS_TOP})
        IF(NOT DEFINED TEST_ENV)
		MESSAGE(FATAL_ERROR "Cannot build with CFI without environment variable FEGS_TOP set which points to the top folder in the CFI installation")
        ENDIF()
	FIND_LIBRARY(CFI_LIBRARY_API NAMES cadfixapi PATHS $ENV{FEGS_LIB})
	FIND_LIBRARY(CFI_LIBRARY_CXX NAMES oocfi_cxx.a PATHS $ENV{FEGS_TOP}/cadfixdev/oocfi/cxx)
        IF(CFI_LIBRARY_API)
		FIND_PATH (CFI_INCLUDE_DIR_HXX cadfixapi.hxx PATHS $ENV{FEGS_TOP}/cadfixdev/oocfi/cxx/cadfixapi)
		FIND_PATH (CFI_INCLUDE_DIR cfiStandardFun.h PATHS $ENV{FEGS_TOP}/cadfixdev/include)

            IF(CFI_INCLUDE_DIR)

		MESSAGE(STATUS "Found CFI Libraries: ${CFI_LIBRARY_API}")

                INCLUDE_DIRECTORIES(NekMeshUtils ${CFI_INCLUDE_DIR_HXX})
                INCLUDE_DIRECTORIES(NekMeshUtils ${CFI_INCLUDE_DIR})

            ELSE()
                MESSAGE(FATAL_ERROR "Cannot find cadfixapi headers")
            ENDIF()
        ELSE()
            MESSAGE(FATAL_ERROR "Cannot find cadfixapi libraries")
        ENDIF()
    ENDIF()
ENDIF()

INCLUDE_DIRECTORIES(SYSTEM ${CFI_INCLUDE_DIR_HXX})
