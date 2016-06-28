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
            TKXSBase TKPShape TKTopAlgo)
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
            URL http://files.opencascade.com/OCCT/OCC_6.9.0_release/opencascade-6.9.0.tgz
            URL_MD5 ba87fe9f5ca47e3dfd62aad7223f0e7f
            STAMP_DIR ${TPBUILD}/stamp
            BINARY_DIR ${TPBUILD}/opencascade-6.9
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/opencascade-6.9
            INSTALL_DIR ${TPBUILD}/opencascade-6.9/dist
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DINSTALL_DIR:PATH=${TPDIST}
                -DCMAKE_CXX_FLAGS:STRING=-DTIXML_USE_STL
                ${TPSRC}/opencascade-6.9
            )

        IF (APPLE)
            SET(OCC_LIBDIR ${TPBUILD}/opencascade-6.9/mac64/clang/)
        ENDIF()
        
        EXTERNALPROJECT_ADD_STEP(opencascade-6.9 post-install
            COMMAND ${CMAKE_COMMAND} -E make_directory ${TPDIST}/lib
            COMMAND ${CMAKE_COMMAND} -E make_directory ${TPDIST}/include
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${TPBUILD}/opencascade-6.9/dist/inc ${TPDIST}/include
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${OCC_LIBDIR} ${TPDIST}/lib
            DEPENDEES install)

        # Patch OS X libraries to fix install name problems.
        #EXTERNALPROJECT_ADD_STEP(opencascade-6.9 patch-install-path
        #    COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/scripts/patch-occ.sh ${TPSRC}/opencascade-6.8/i686/lib ${CMAKE_INSTALL_PREFIX}/${NEKTAR_LIB_DIR}
        #    DEPENDEES build
        #    DEPENDERS install)

        MESSAGE(STATUS "Build OpenCascade: ${TPDIST}/lib")
        LINK_DIRECTORIES(${TPDIST}/lib)
        INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include)
    ELSE()
        ADD_CUSTOM_TARGET(opencascade-6.9 ALL)
        MESSAGE(STATUS "Found OpenCASCADE: ${OCC_LIBRARY_DIR}")
        SET(OPENCASCADE_CONFIG_INCLUDE_DIR ${OCC_INCLUDE_DIR})
    ENDIF()
ENDIF()
