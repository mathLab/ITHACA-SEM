########################################################################
#
# ThirdParty configuration for Nektar++
#
# LMFIT
#
########################################################################

IF(NEKTAR_USE_MESH)


SET(BUILD_LMFIT ON)

OPTION(THIRDPARTY_BUILD_LMFIT
    "Build lmfit library from ThirdParty." ${BUILD_LMFIT})

IF (THIRDPARTY_BUILD_LMFIT)
    INCLUDE(ExternalProject)
    EXTERNALPROJECT_ADD(
        lmfit-5.1
        PREFIX ${TPSRC}
        URL http://ae-nektar.ae.ic.ac.uk/~mt4313/lmfit.zip
        URL_MD5 7a32f672fbe911ecc8302eef385b3780
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/lmfit-5.1
        BINARY_DIR ${TPBUILD}/lmfit-5.1
        TMP_DIR ${TPBUILD}/lmfit-5.1-tmp
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            ${TPSRC}/lmfit-5.1
    )
    SET(LMFIT_LIBRARY lmfit CACHE FILEPATH
        "lmfit library" FORCE)
    SET(LMFIT_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
        "lmfit include" FORCE)

    LINK_DIRECTORIES(${TPDIST}/lib)

    IF (WIN32)
        MESSAGE(STATUS
                "Build lmfit: ${TPDIST}/${LIB_DIR}/${LMFIT_LIBRARY}.dll")
    ELSE ()
        MESSAGE(STATUS
                "Build lmfit: ${TPDIST}/${LIB_DIR}/lib${LMFIT_LIBRARY}.a")
    ENDIF ()

    SET(LMFIT_CONFIG_INCLUDE_DIR ${TPINC})
ELSE()
    ADD_CUSTOM_TARGET(lmfit-5.1 ALL)
    MESSAGE(STATUS "Found lmfit: ${LMFIT_LIBRARY}")
    SET(LMFIT_CONFIG_INCLUDE_DIR ${LMFIT_INCLUDE_DIR})
ENDIF (THIRDPARTY_BUILD_LMFIT)

INCLUDE_DIRECTORIES(SYSTEM ${LMFIT_INCLUDE_DIR})

ENDIF(NEKTAR_USE_MESH)
