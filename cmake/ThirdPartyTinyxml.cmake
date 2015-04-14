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
FIND_PATH   (TINYXML_INCLUDE_DIR tinyxml.h PATHS /opt/local/include)
FIND_LIBRARY(TINYXML_LIBRARY NAMES "tinyxml" PATHS /opt/local/lib)

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
    SET(TINYXML_LIBRARY tinyxml CACHE FILEPATH
        "TinyXML library" FORCE)
    SET(TINYXML_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
        "TinyXML include" FORCE)

    LINK_DIRECTORIES(${TPDIST}/lib)

    IF (WIN32)
        MESSAGE(STATUS 
                "Build TinyXML: ${TPDIST}/${LIB_DIR}/${TINYXML_LIBRARY}.dll")
    ELSE ()
        MESSAGE(STATUS 
                "Build TinyXML: ${TPDIST}/${LIB_DIR}/lib${TINYXML_LIBRARY}.a")
    ENDIF ()

    SET(TINYXML_CONFIG_INCLUDE_DIR ${TPINC})
ELSE()
    ADD_CUSTOM_TARGET(tinyxml-2.6.2 ALL)
    MESSAGE(STATUS "Found TinyXML: ${TINYXML_LIBRARY}")
    SET(TINYXML_CONFIG_INCLUDE_DIR ${TINYXML_INCLUDE_DIR})
ENDIF (THIRDPARTY_BUILD_TINYXML)

INCLUDE_DIRECTORIES(SYSTEM ${TINYXML_INCLUDE_DIR})

MARK_AS_ADVANCED(TINYXML_INCLUDE_DIR)
MARK_AS_ADVANCED(TINYXML_LIBRARY)
