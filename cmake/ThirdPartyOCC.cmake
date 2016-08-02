########################################################################
#
# ThirdParty configuration for Nektar++
#
# OpenCascade
#
########################################################################

IF(NEKTAR_USE_MESHGEN)
    # Try to find installed version of OpenCascade
    INCLUDE(FindOCC)

    IF (OCC_FOUND)
        SET(BUILD_OCC OFF)
    ELSE()
        SET(BUILD_OCC ON)
    ENDIF()
    
    OPTION(THIRDPARTY_BUILD_OCC "Build OpenCascade library from ThirdParty."
        ${BUILD_OCC})

    IF (THIRDPARTY_BUILD_OCC)
        INCLUDE(ExternalProject)

        SET(OCC_LIBRARIES_TMP PTKernel TKernel TKMath TKBRep TKIGES TKSTEP TKSTEPAttr
            TKSTEP209 TKSTEPBase TKShapeSchema TKGeomBase TKGeomAlgo TKG3d TKG2d
            TKXSBase TKPShape TKTopAlgo TKShHealing)
        FOREACH(OCC_LIB ${OCC_LIBRARIES_TMP})
            LIST(APPEND OCC_LIBRARIES ${TPDIST}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${OCC_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
        ENDFOREACH()
        UNSET(OCC_LIBRARIES_TMP)

        IF(WIN32)
            MESSAGE(SEND_ERROR "Cannot currently use OpenCascade with Nektar++ on Windows")
        ENDIF()
        
        EXTERNALPROJECT_ADD(
            opencascade-6.9
            PREFIX ${TPSRC}
            URL http://ae-nektar.ae.ic.ac.uk/~dmoxey/OCE-0.17.2.tar.gz
            URL_MD5 bf2226be4cd192606af677cf178088e5
            STAMP_DIR ${TPBUILD}/stamp
            BINARY_DIR ${TPBUILD}/opencascade-6.9
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/opencascade-6.9
            INSTALL_DIR ${TPBUILD}/opencascade-6.9/dist
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DOCE_INSTALL_PREFIX:PATH=${TPDIST}
                -DOCE_TESTING=OFF
                -DOCE_VISUALISATION=OFF
                -DOCE_DISABLE_X11=ON
                ${TPSRC}/opencascade-6.9
            )

        # Patch OS X libraries to fix install name problems.
        EXTERNALPROJECT_ADD_STEP(opencascade-6.9 patch-install-path
            COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/scripts/patch-occ.sh ${TPDIST}/lib ${CMAKE_INSTALL_PREFIX}/${NEKTAR_LIB_DIR}
            ALWAYS 1
            DEPENDEES install)

        MESSAGE(STATUS "Build OpenCascade: ${TPDIST}/lib")
        LINK_DIRECTORIES(${TPDIST}/lib)
        INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include/oce)
    ELSE()
        ADD_CUSTOM_TARGET(opencascade-6.9 ALL)
        SET(OPENCASCADE_CONFIG_INCLUDE_DIR ${OCC_INCLUDE_DIR})
        INCLUDE_DIRECTORIES(SYSTEM ${OCC_INCLUDE_DIR})
    ENDIF()
ENDIF()
