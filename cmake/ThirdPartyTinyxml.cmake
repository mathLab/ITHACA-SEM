########################################################################
#
# ThirdParty configuration for Nektar++
#
# TinyXML
#
########################################################################

OPTION(NEKTAR_USE_TINYXML_STL "Use STL with TinyXML library." ON)
MARK_AS_ADVANCED(NEKTAR_USE_TINYXML_STL)

# First search for system TinyXML installs. Hint /opt/local for MacPorts.
FIND_PATH   (TINYXML_INCLUDE_DIR tinyxml.h)
FIND_LIBRARY(TINYXML_LIBRARY NAMES "tinyxml")

# If we have our library then don't build TinyXML.
IF (TINYXML_INCLUDE_DIR AND TINYXML_LIBRARY)
    SET(BUILD_TINYXML OFF)
ELSE()
    SET(BUILD_TINYXML ON)
ENDIF ()

OPTION(THIRDPARTY_BUILD_TINYXML
    "Build TinyXML library from ThirdParty." ${BUILD_TINYXML})

IF (THIRDPARTY_BUILD_TINYXML)
    INCLUDE(ExternalProject)

    find_program(HAS_PATCH patch)

    IF(HAS_PATCH)
        EXTERNALPROJECT_ADD(
            tinyxml-2.6.2
            PREFIX ${TPSRC}
            URL ${TPURL}/tinyxml_2_6_2.tar.bz2
            URL_MD5 240beaeb45f63b154c9801eef7561eac
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/tinyxml-2.6.2
            BINARY_DIR ${TPBUILD}/tinyxml-2.6.2
            TMP_DIR ${TPBUILD}/tinyxml-2.6.2-tmp
            INSTALL_DIR ${TPDIST}
            PATCH_COMMAND patch -d ${TPSRC}/tinyxml-2.6.2 -o tmp < ${CMAKE_SOURCE_DIR}/cmake/scripts/tinyxml.patch
            COMMAND ${CMAKE_COMMAND} -E copy ${TPSRC}/tinyxml-2.6.2/tmp ${TPSRC}/tinyxml-2.6.2/tinyxmlparser.cpp
            COMMAND ${CMAKE_COMMAND} -E remove ${TPSRC}/tinyxml-2.6.2/tmp
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                -DCMAKE_CXX_FLAGS:STRING=-DTIXML_USE_STL
                ${TPSRC}/tinyxml-2.6.2
            )
    ELSE()
        MESSAGE(STATUS "patch utility not found, cannot apply patch to tinyxml, Nektar++ will still function")
        EXTERNALPROJECT_ADD(
            tinyxml-2.6.2
            PREFIX ${TPSRC}
            URL ${TPURL}/tinyxml_2_6_2.tar.bz2
            URL_MD5 240beaeb45f63b154c9801eef7561eac
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/tinyxml-2.6.2
            BINARY_DIR ${TPBUILD}/tinyxml-2.6.2
            TMP_DIR ${TPBUILD}/tinyxml-2.6.2-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                -DCMAKE_CXX_FLAGS:STRING=-DTIXML_USE_STL
                ${TPSRC}/tinyxml-2.6.2
            )
    ENDIF()

    THIRDPARTY_LIBRARY(TINYXML_LIBRARY STATIC tinyxml DESCRIPTION "TinyXML library")
    SET(TINYXML_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
        "TinyXML include" FORCE)
    MESSAGE(STATUS "Build TinyXML: ${TINYXML_LIBRARY}")
    SET(TINYXML_CONFIG_INCLUDE_DIR ${TPINC})
ELSE()
    ADD_CUSTOM_TARGET(tinyxml-2.6.2 ALL)
    MESSAGE(STATUS "Found TinyXML: ${TINYXML_LIBRARY}")
    SET(TINYXML_CONFIG_INCLUDE_DIR ${TINYXML_INCLUDE_DIR})
ENDIF (THIRDPARTY_BUILD_TINYXML)

INCLUDE_DIRECTORIES(${TINYXML_INCLUDE_DIR})

MARK_AS_ADVANCED(TINYXML_INCLUDE_DIR)
MARK_AS_ADVANCED(TINYXML_LIBRARY)
MARK_AS_ADVANCED(TINYXML_CONFIG_INCLUDE_DIR)
