########################################################################
#
# ThirdParty configuration for Nektar++
#
# Levmar
#
########################################################################



IF(NEKTAR_USE_MESH)


SET(BUILD_LEVMAR ON)


OPTION(THIRDPARTY_BUILD_LEVMAR
    "Build Levmar library from ThirdParty." ${BUILD_LEVMAR})

IF (THIRDPARTY_BUILD_LEVMAR)
    INCLUDE(ExternalProject)
    EXTERNALPROJECT_ADD(
        levmar-5.2
        PREFIX ${TPSRC}
        URL http://ae-nektar.ae.ic.ac.uk/~mt4313/levmar.zip
        URL_MD5 7139c9790e3ed4cb5fe2d5be6b1e30d5
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/levmar-5.2
        BINARY_DIR ${TPBUILD}/levmar-5.2
        TMP_DIR ${TPBUILD}/levmar-5.2-tmp
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            ${TPSRC}/levmar-5.2
    )
    SET(LEVMAR_LIBRARY levmar CACHE FILEPATH
        "Levmar library" FORCE)
    SET(LEVMAR_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
        "Levmar include" FORCE)

    LINK_DIRECTORIES(${TPDIST}/lib)

    IF (WIN32)
        MESSAGE(STATUS
                "Build Levmar: ${TPDIST}/${LIB_DIR}/${LEVMAR_LIBRARY}.dll")
    ELSE ()
        MESSAGE(STATUS
                "Build Levmar: ${TPDIST}/${LIB_DIR}/lib${LEVMAR_LIBRARY}.a")
    ENDIF ()

    SET(TRIANGLE_CONFIG_INCLUDE_DIR ${TPINC})
ELSE()
    ADD_CUSTOM_TARGET(levmar-5.2 ALL)
    MESSAGE(STATUS "Found Levmar: ${LEVMAR_LIBRARY}")
    SET(TRIANGLE_CONFIG_INCLUDE_DIR ${LEVMAR_INCLUDE_DIR})
ENDIF (THIRDPARTY_BUILD_LEVMAR)

INCLUDE_DIRECTORIES(SYSTEM ${LEVMAR_INCLUDE_DIR})

ENDIF(NEKTAR_USE_MESH)
