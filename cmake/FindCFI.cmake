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

        SET(TEST_ENV $ENV{CFI_DIR})
        IF(NOT DEFINED TEST_ENV)
            MESSAGE(FATAL_ERROR "Cannot build with CFI without environment variable CFI_DIR set which points to cadfix1100fcs folder in the CFI installation")
        ENDIF()

        FIND_LIBRARY(CFI_LIBRARY_API NAMES cadfixapi PATHS $ENV{CFI_DIR}/lib64)

        IF(CFI_LIBRARY_API)
            FIND_PATH (CFI_INCLUDE_DIR_HXX cadfixapi.hxx PATHS $ENV{CFI_DIR}/oocfi/cxx/cadfixapi)
            FIND_PATH (CFI_INCLUDE_DIR cfiStandardFun.h PATHS $ENV{CFI_DIR}/include)

            IF(CFI_INCLUDE_DIR)
                SET(CFI_LIBRARIES_TMP cadfixapi extra)
                FOREACH(CFI_LIBRARIES_TMP ${CFI_LIBRARIES_TMP})
                    LIST(APPEND CFI_LIBRARIES $ENV{CFI_DIR}/lib64/lib${CFI_LIBRARIES_TMP}.so)
                ENDFOREACH()

                MESSAGE(STATUS "cfi libraries: ${CFI_LIBRARIES}")

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
