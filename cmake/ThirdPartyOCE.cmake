########################################################################
#
# ThirdParty configuration for Nektar++
#
# OpenCascade
#
########################################################################

IF(NEKTAR_USE_MESHGEN)
    # Try to find installed version of OpenCascade
    INCLUDE(FindOCCE)

    IF (OCE_FOUND)
        SET(BUILD_OCE OFF)
    ELSE()
        SET(BUILD_OCE ON)
    ENDIF()

    OPTION(THIRDPARTY_BUILD_OCE "Build OpenCascade community edition library from ThirdParty."
        ${BUILD_OCE})

    IF (THIRDPARTY_BUILD_OCE)
        INCLUDE(ExternalProject)

        SET(OCE_LIBRARIES_TMP PTKernel TKernel TKMath TKBRep TKIGES TKSTEP TKSTEPAttr
            TKSTEP209 TKSTEPBase TKShapeSchema TKGeomBase TKGeomAlgo TKG3d TKG2d
            TKXSBase TKPShape TKTopAlgo TKShHealing)
        FOREACH(OCE_LIB ${OCE_LIBRARIES_TMP})
            LIST(APPEND OCE_LIBRARIES ${TPDIST}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${OCE_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
        ENDFOREACH()
        UNSET(OCE_LIBRARIES_TMP)

        IF(WIN32)
            MESSAGE(SEND_ERROR "Cannot currently use OpenCascade with Nektar++ on Windows")
        ENDIF()

        EXTERNALPROJECT_ADD(
            oce-0.17
            PREFIX ${TPSRC}
            URL ${TPURL}/OCE-0.17.2.tar.gz
            URL_MD5 bf2226be4cd192606af677cf178088e5
            STAMP_DIR ${TPBUILD}/stamp
            BINARY_DIR ${TPBUILD}/oce-0.17
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/oce-0.17
            INSTALL_DIR ${TPBUILD}/oce-0.17/dist
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DOCE_INSTALL_PREFIX:PATH=${TPDIST}
                -DOCE_TESTING=OFF
                -DOCE_VISUALISATION=OFF
                -DOCE_DISABLE_X11=ON
                -DOCE_OCAF=OFF
                ${TPSRC}/oce-0.17
            )

        # Patch OS X libraries to fix install name problems.
        EXTERNALPROJECT_ADD_STEP(oce-0.17 patch-install-path
            COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/scripts/patch-occ.sh ${TPDIST}/lib ${CMAKE_INSTALL_PREFIX}/${NEKTAR_LIB_DIR}
            ALWAYS 1
            DEPENDEES install)

        MESSAGE(STATUS "Build OpenCascade community edition: ${TPDIST}/lib")
        LINK_DIRECTORIES(${TPDIST}/lib)
        INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include/oce)
    ELSE()
        ADD_CUSTOM_TARGET(oce-0.17 ALL)
        SET(OPENCASCADE_CONFIG_INCLUDE_DIR ${OCE_INCLUDE_DIRS})
        INCLUDE_DIRECTORIES(SYSTEM ${OCE_INCLUDE_DIRS})
    ENDIF()
ENDIF()
