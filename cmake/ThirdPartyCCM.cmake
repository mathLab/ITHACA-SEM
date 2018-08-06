########################################################################
#
# ThirdParty configuration for Nektar++
#
# Star CCM i/o
#
########################################################################
OPTION(NEKTAR_USE_CCM
   "use CCM star i/o" OFF)

IF(NEKTAR_USE_CCM)

# First search for system ccmioL installs. Hint /usr/local 
FIND_PATH   (CCMIO_INCLUDE_DIR ccmio.h PATHS /usr/local/include ${CCM_DIR} PATH_SUFFIXES libccmio)
FIND_LIBRARY(CCMIO_LIBRARY NAMES "ccmio" PATHS /usr/local/lib ${CCM_DIR} PATH_SUFFIXES lib)

# If we have our library then don't build CCMIO.
IF (CCMIO_INCLUDE_DIR AND CCMIO_LIBRARY)
    SET(BUILD_CCMIO OFF)
ELSE()
    SET(BUILD_CCMIO ON)
ENDIF ()

OPTION(THIRDPARTY_BUILD_CCMIO 
    "Build CCMIO library from ThirdParty if permitted." ${BUILD_CCMIO})

IF (THIRDPARTY_BUILD_CCMIO)
    INCLUDE(ExternalProject)
    MESSAGE(WARNING "We are seeking permission to distribute ccmio with Nektar++. If you are entitled to use libccmio please contact nektar-users@imperial.ac.uk and place the file ccmio-2.06.tar.bz2 in the director $NEKTAR/ThirdParty")
    EXTERNALPROJECT_ADD(
        ccmio-2.06
        PREFIX ${TPSRC}
        URL ${TPURL}/ccmio-2.06.tar.bz2
        URL_MD5 809ee34a983cbc8931ca23879d92b4d0
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/ccmio-2.06
        BINARY_DIR ${TPBUILD}/ccmio-2.06
        TMP_DIR ${TPBUILD}/ccmio-2.06-tmp
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            ${TPSRC}/ccmio-2.06
    )
    SET(CCMIO_LIBRARY ccmio CACHE FILEPATH
        "CCMIO library" FORCE)
    SET(CCMIO_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
        "CCMIO include" FORCE)

    LINK_DIRECTORIES(${TPDIST}/lib)

    INCLUDE_DIRECTORIES(NekMesh ${CCMIO_INCLUDE_DIR})

    IF (WIN32)
        MESSAGE(STATUS 
                "Build CCMIO: ${TPDIST}/${LIB_DIR}/${CCMIO_LIBRARY}.dll")
    ELSE ()
        MESSAGE(STATUS 
                "Build CCMIO: ${TPDIST}/${LIB_DIR}/lib${CCMIO_LIBRARY}.a")
    ENDIF ()

    SET(CCMIO_CONFIG_INCLUDE_DIR ${TPINC})

    set(CCMIO_LIBRARIES
        ccmio
        adf
    )
ELSE()
    ADD_CUSTOM_TARGET(ccmio-2.06 ALL)
    MESSAGE(STATUS "Found CCMIO: ${CCMIO_LIBRARY}")
    SET(CCMIO_CONFIG_INCLUDE_DIR ${CCMIO_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(NekMesh ${CCMIO_INCLUDE_DIR})
    LINK_DIRECTORIES(${CCMIO_LIBRARY_DIR})

ENDIF (THIRDPARTY_BUILD_CCMIO)

SET(CCMIO_LIBRARIES ccmio adf)

MARK_AS_ADVANCED(CCMIO_INCLUDE_DIR)
MARK_AS_ADVANCED(CCMIO_LIBRARY)

ENDIF(NEKTAR_USE_CCM)
