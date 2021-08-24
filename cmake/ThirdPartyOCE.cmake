########################################################################
#
# ThirdParty configuration for Nektar++
#
# OpenCascade
#
########################################################################

IF(NEKTAR_USE_MESHGEN)
    #required opencascade libraries
    SET(OCC_LIB_LIST
        TKFillet
        TKMesh
        TKernel
        TKG2d
        TKG3d
        TKMath
        TKIGES
        TKSTL
        TKShHealing
        TKXSBase
        TKBool
        TKBO
        TKBRep
        TKTopAlgo
        TKGeomAlgo
        TKGeomBase
        TKOffset
        TKPrim
        TKSTEP
        TKSTEPBase
        TKSTEPAttr
        TKHLR
        TKFeat
        TKXCAF
        TKLCAF
        TKXDESTEP
    )

    # Try to find installed version of OpenCascade
    INCLUDE(FindOCC)

    IF (OCC_FOUND)
        SET(BUILD_OCE OFF)
    ELSE()
        SET(BUILD_OCE ON)
    ENDIF()

    OPTION(THIRDPARTY_BUILD_OCE "Build OpenCascade community edition library from ThirdParty."
        ${BUILD_OCE})

    IF (THIRDPARTY_BUILD_OCE)
        INCLUDE(ExternalProject)

        IF(WIN32)
            MESSAGE(SEND_ERROR "Cannot currently use OpenCascade with Nektar++ on Windows")
        ENDIF()

        EXTERNALPROJECT_ADD(
            oce-0.18.3
            PREFIX ${TPSRC}
            URL ${TPURL}/OCE-0.18.3.tar.gz
            URL_MD5 1686393c8493bbbb2f3f242330b33cba
            STAMP_DIR ${TPBUILD}/stamp
            BINARY_DIR ${TPBUILD}/oce-0.18.3
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/oce-0.18.3
            INSTALL_DIR ${TPBUILD}/oce-0.18.3/dist
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DOCE_INSTALL_PREFIX:PATH=${TPDIST}
                -DOCE_TESTING=OFF
                -DOCE_VISUALISATION=OFF
                ${TPSRC}/oce-0.18.3
            )

        # Patch OS X libraries to fix install name problems.
        IF(APPLE)
            EXTERNALPROJECT_ADD_STEP(oce-0.18.3 patch-install-path
                COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/scripts/patch-occ.sh ${TPDIST}/lib ${CMAKE_INSTALL_PREFIX}/${NEKTAR_LIB_DIR}
                DEPENDEES install)
        ENDIF()

        THIRDPARTY_LIBRARY(OCC_LIBRARIES SHARED ${OCC_LIB_LIST} DESCRIPTION "OpenCascade libs")
        SET(OCC_INCLUDE_DIR ${TPDIST}/include/oce CACHE FILEPATH "OCC include" FORCE)
        MESSAGE(STATUS "Build OpenCascade community edition: ${TPDIST}/lib")
    ELSE()
        ADD_CUSTOM_TARGET(oce-0.18.3 ALL)
    ENDIF()
ENDIF()

INCLUDE_DIRECTORIES(SYSTEM ${OCC_INCLUDE_DIR})
