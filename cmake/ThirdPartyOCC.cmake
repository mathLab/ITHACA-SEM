########################################################################
#
# ThirdParty configuration for Nektar++
#
# OpenCascade
#
########################################################################

IF(NEKTAR_USE_MESHGEN)
    SET(BUILD_OCC ON)
    OPTION(THIRDPARTY_DOWNLOAD_OCC
        "Download pre-compiled versions of OpenCascade." ${BUILD_OCC})

    IF (THIRDPARTY_DOWNLOAD_OCC)
        INCLUDE(ExternalProject)

        SET(OCC_LIBS_TMP PTKernel TKernel TKMath TKBRep TKIGES TKSTEP TKSTEPAttr
            TKSTEP209 TKSTEPBase TKShapeSchema TKGeomBase TKGeomAlgo TKG3d TKG2d
            TKXSBase TKPShape TKTopAlgo)
        FOREACH(OCC_LIB ${OCC_LIBS_TMP})
            LIST(APPEND OCC_LIBS ${TPDIST}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${OCC_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
        ENDFOREACH()

        IF(WIN32)
            MESSAGE(SEND_ERROR "Cannot currently use OpenCascade with Nektar++ on Windows")
        ELSEIF(APPLE)
            EXTERNALPROJECT_ADD(
                opencascade-6.8
                PREFIX ${TPSRC}
                URL ${TPURL}/OCC680osx64.tgz
                URL_MD5 626292523b0691304f0fa271989fbc44
                STAMP_DIR ${TPBUILD}/stamp
                BINARY_DIR ${TPBUILD}/opencascade-6.8
                DOWNLOAD_DIR ${TPSRC}
                SOURCE_DIR ${TPSRC}/opencascade-6.8
                INSTALL_DIR ${TPDIST}
                UPDATE_COMMAND ""
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                INSTALL_COMMAND cp -a ${TPSRC}/opencascade-6.8/i686/lib/. ${TPDIST}/lib/ COMMAND cp -a ${TPSRC}/opencascade-6.8/i686/inc/. ${TPDIST}/include/
                )

            # Patch OS X libraries to fix install name problems.
            EXTERNALPROJECT_ADD_STEP(opencascade-6.8 patch-install-path
                COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/scripts/patch-occ.sh ${TPSRC}/opencascade-6.8/i686/lib ${CMAKE_INSTALL_PREFIX}/${NEKTAR_LIB_DIR}
                DEPENDEES build
                DEPENDERS install)
        ELSE()
            EXTERNALPROJECT_ADD(
                opencascade-6.8
                PREFIX ${TPSRC}
                URL ${TPURL}/OCC680lin64.tgz
                URL_MD5 d655b6f50998bb9600e081907c247793
                STAMP_DIR ${TPBUILD}/stamp
                DOWNLOAD_DIR ${TPSRC}
                SOURCE_DIR ${TPSRC}/opencascade-6.8
                INSTALL_DIR ${TPDIST}
                UPDATE_COMMAND ""
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                INSTALL_COMMAND cp -a ${TPSRC}/opencascade-6.8/lib/. ${TPDIST}/lib/ COMMAND cp -a ${TPSRC}/opencascade-6.8/inc/. ${TPDIST}/include/
                )
        ENDIF()

        MESSAGE(STATUS "Build OpenCascade: download binaries")
        LINK_DIRECTORIES(${TPDIST}/lib)
        INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include)
    ENDIF (THIRDPARTY_DOWNLOAD_OCC)
ENDIF(NEKTAR_USE_MESHGEN)
