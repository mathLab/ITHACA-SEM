########################################################################
#
# ThirdParty configuration for Nektar++
#
# TRIANGLE
#
########################################################################

IF(NEKTAR_USE_MESHGEN)
    # Search for system-installed Triangle installation
    FIND_LIBRARY(TRIANGLE_LIBRARY NAMES triangle)
    FIND_PATH(TRIANGLE_INCLUDE_DIR triangle.h)

    IF(TRIANGLE_LIBRARY AND TRIANGLE_INCLUDE_DIR)
        SET(BUILD_TRIANGLE OFF)
    ELSE()
        SET(BUILD_TRIANGLE ON)
    ENDIF()

    OPTION(THIRDPARTY_BUILD_TRIANGLE
    "Build Triangle library from ThirdParty." ${BUILD_TRIANGLE})

    IF (THIRDPARTY_BUILD_TRIANGLE)
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            triangle-1.6
            PREFIX ${TPSRC}
            URL ${TPURL}/triangle.zip
            URL_MD5 357cb7107f51f3f89940c47435d4fa49
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/triangle-1.6
            BINARY_DIR ${TPBUILD}/triangle-1.6
            TMP_DIR ${TPBUILD}/triangle-1.6-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                ${TPSRC}/triangle-1.6
            )
        THIRDPARTY_LIBRARY(TRIANGLE_LIBRARY STATIC triangle
            DESCRIPTION "Triangle library")
        SET(TRIANGLE_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "Triangle include" FORCE)
        MESSAGE(STATUS "Build Triangle: ${TRIANGLE_LIBRARY}")
        SET(TRIANGLE_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE()
        ADD_CUSTOM_TARGET(triangle-1.6 ALL)
        MESSAGE(STATUS "Found Triangle: ${TRIANGLE_LIBRARY}")
        SET(TRIANGLE_CONFIG_INCLUDE_DIR ${TRIANGLE_INCLUDE_DIR})
    ENDIF (THIRDPARTY_BUILD_TRIANGLE)

    MARK_AS_ADVANCED(TRIANGLE_LIBRARY)
    MARK_AS_ADVANCED(TRIANGLE_INCLUDE_DIR)
    INCLUDE_DIRECTORIES(${TRIANGLE_INCLUDE_DIR})
ENDIF(NEKTAR_USE_MESHGEN)
