########################################################################
#
# ThirdParty configuration for Nektar++
#
# TETGEN
#
########################################################################

IF(NEKTAR_USE_MESHGEN)
    # Search for system-installed Triangle installation
    FIND_LIBRARY(TETGEN_LIBRARY NAMES tet)
    FIND_PATH(TETGEN_INCLUDE_DIR tetgen.h)

    IF(TETGEN_LIBRARY AND TETGEN_INCLUDE_DIR)
        SET(BUILD_TETGEN OFF)
    ELSE()
        SET(BUILD_TETGEN ON)
    ENDIF()

    OPTION(THIRDPARTY_BUILD_TETGEN
        "Build TetGen library from ThirdParty." ${BUILD_TETGEN})

    IF (THIRDPARTY_BUILD_TETGEN)
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            tetgen-1.5
            PREFIX ${TPSRC}
            URL ${TPURL}/tetgen.zip
            URL_MD5 6d62e63f9b1e7a8ce53d5bc87e6a0a09
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/tetgen-1.5
            BINARY_DIR ${TPBUILD}/tetgen-1.5
            TMP_DIR ${TPBUILD}/tetgen-1.5-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            ${TPSRC}/tetgen-1.5
            )
        THIRDPARTY_LIBRARY(TETGEN_LIBRARY STATIC tetgen
            DESCRIPTION "Tetgen library")
        SET(TETGEN_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "TetGen include" FORCE)
        MESSAGE(STATUS "Build TetGen: ${TETGEN_LIBRARY}")
        SET(TETGEN_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE()
        ADD_CUSTOM_TARGET(tetgen-1.5 ALL)
        MESSAGE(STATUS "Found TetGen: ${TETGEN_LIBRARY}")
        SET(TETGEN_CONFIG_INCLUDE_DIR ${TETGEN_INCLUDE_DIR})
    ENDIF()

    MARK_AS_ADVANCED(TETGEN_LIBRARY)
    MARK_AS_ADVANCED(TETGEN_INCLUDE_DIR)
    INCLUDE_DIRECTORIES(${TETGEN_INCLUDE_DIR})
ENDIF(NEKTAR_USE_MESHGEN)
