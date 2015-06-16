########################################################################
#
# ThirdParty configuration for Nektar++
#
# Triangle
#
########################################################################

# First search for system Triangle installs. Hint /opt/local for MacPorts.
FIND_PATH   (TRIANGLE_INCLUDE_DIR triangle.h PATHS /opt/local/include)
FIND_LIBRARY(TRIANGLE_LIBRARY NAMES "triangle" PATHS /opt/local/lib)

# If we have our library then don't build Triangle
IF (TRIANGLE_INCLUDE_DIR AND TRIANGLE_LIBRARY)
    SET(BUILD_TRIANGLE OFF)
ELSE()
    SET(BUILD_TRIANGLE ON)
ENDIF ()

OPTION(THIRDPARTY_BUILD_TRIANGLE
    "Build Triangle library from ThirdParty." ${BUILD_TRIANGLE})

IF (THIRDPARTY_BUILD_TRIANGLE)
    INCLUDE(ExternalProject)
    EXTERNALPROJECT_ADD(
        triangle-1.6
        PREFIX ${TPSRC}
        URL ${TPURL}/triangle.zip
        URL_MD5 ebc48aa34d79b2a03080c8647e8164a6
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
    SET(TRIANGLE_LIBRARY triangle CACHE FILEPATH
        "Triangle library" FORCE)
    SET(TRIANGLE_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
        "Triangle include" FORCE)

    LINK_DIRECTORIES(${TPDIST}/lib)

    IF (WIN32)
        MESSAGE(STATUS 
                "Build Triangle: ${TPDIST}/${LIB_DIR}/${TRIANGLE_LIBRARY}.dll")
    ELSE ()
        MESSAGE(STATUS 
                "Build Triangle: ${TPDIST}/${LIB_DIR}/lib${TRIANGLE_LIBRARY}.a")
    ENDIF ()

    SET(TRIANGLE_CONFIG_INCLUDE_DIR ${TPINC})
ELSE()
    ADD_CUSTOM_TARGET(triangle-1.6 ALL)
    MESSAGE(STATUS "Found Triangle: ${TRIANGLE_LIBRARY}")
    SET(TRIANGLE_CONFIG_INCLUDE_DIR ${TRIANGLE_INCLUDE_DIR})
ENDIF (THIRDPARTY_BUILD_TRIANGLE)

INCLUDE_DIRECTORIES(SYSTEM ${TRIANGLE_INCLUDE_DIR})

