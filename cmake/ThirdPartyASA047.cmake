########################################################################
#
# ThirdParty configuration for Nektar++
#
# ASA047
#
########################################################################



IF(NEKTAR_USE_MESH)


SET(BUILD_ASA ON)


OPTION(THIRDPARTY_BUILD_ASA
    "Build Asa047 library from ThirdParty." ${BUILD_ASA})

IF (THIRDPARTY_BUILD_ASA)
    INCLUDE(ExternalProject)
    EXTERNALPROJECT_ADD(
        asa047-1.0
        PREFIX ${TPSRC}
        URL http://ae-nektar.ae.ic.ac.uk/~mt4313/asa.zip
        URL_MD5 e5b120d1b8759cc265788a12aef0b0b9
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/asa047-1.0
        BINARY_DIR ${TPBUILD}/asa047-1.0
        TMP_DIR ${TPBUILD}/asa047-1.0-tmp
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            ${TPSRC}/asa047-1.0
    )
    SET(ASA_LIBRARY asa047 CACHE FILEPATH
        "Asa047 library" FORCE)
    SET(ASA_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
        "Asa047 include" FORCE)

    LINK_DIRECTORIES(${TPDIST}/lib)

    IF (WIN32)
        MESSAGE(STATUS
                "Build Asa047: ${TPDIST}/${LIB_DIR}/${ASA_LIBRARY}.dll")
    ELSE ()
        MESSAGE(STATUS
                "Build Asa047: ${TPDIST}/${LIB_DIR}/lib${ASA_LIBRARY}.a")
    ENDIF ()

    SET(ASA_CONFIG_INCLUDE_DIR ${TPINC})
ELSE()
    ADD_CUSTOM_TARGET(asa047-1.0 ALL)
    MESSAGE(STATUS "Found asa: ${ASA_LIBRARY}")
    SET(TRIANGLE_CONFIG_INCLUDE_DIR ${ASA_INCLUDE_DIR})
ENDIF (THIRDPARTY_BUILD_ASA)

INCLUDE_DIRECTORIES(SYSTEM ${ASA_INCLUDE_DIR})

ENDIF(NEKTAR_USE_MESH)
