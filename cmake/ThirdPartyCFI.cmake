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

        FIND_LIBRARY(CFI_LIBRARY_API NAMES cadfixapi PATHS /home/mike/CFI11/cadfix1100fcs/lib64)

        IF(CFI_LIBRARY_API)
            FIND_PATH (CFI_INCLUDE_DIR_HXX cadfixapi.hxx PATHS /home/mike/CFI11/cadfix1100fcs/oocfi/cxx/cadfixapi)
            FIND_PATH (CFI_INCLUDE_DIR cfiStandardFun.h PATHS /home/mike/CFI11/cadfix1100fcs/include)

            IF(CFI_INCLUDE_DIR)
                GET_FILENAME_COMPONENT(CFI_LIBRARY_PATH ${CFI_LIBRARY_API} PATH)

                SET(CFI_LIBRARIES_TMP cadfixapi extra)
                FOREACH(CFI_LIBRARIES_TMP ${CFI_LIBRARIES_TMP})
                    LIST(APPEND CFI_LIBRARIES ${CFI_LIBRARY_PATH}/lib${CFI_LIBRARIES_TMP}.so)
                ENDFOREACH()

                MESSAGE(STATUS "cfi libraries: ${CFI_LIBRARIES}")
                MESSAGE(STATUS "cfi include: ${CFI_INCLUDE_DIR}")

                LINK_DIRECTORIES(${CFI_LIBRARY_PATH})
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
