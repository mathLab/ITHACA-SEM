########################################################################
#
# ThirdParty configuration for Nektar++
#
# Modified METIS library
#
########################################################################

INCLUDE(ExternalProject)

EXTERNALPROJECT_ADD(
    modmetis-5.1.0
    PREFIX ${TPSRC}
    URL ${TPURL}/modmetis-5.1.0_2.tar.bz2
    URL_MD5 "8a1f1afd39b46a4477c1ea15464cdf89"
    STAMP_DIR ${TPBUILD}/stamp
    DOWNLOAD_DIR ${TPSRC}
    SOURCE_DIR ${TPSRC}/modmetis-5.1.0
    BINARY_DIR ${TPBUILD}/modmetis-5.1.0
    TMP_DIR ${TPBUILD}/modmetis-5.1.0-tmp
    INSTALL_DIR ${TPDIST}
    CONFIGURE_COMMAND ${CMAKE_COMMAND}
        -G ${CMAKE_GENERATOR}
        -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
        -DCMAKE_C_FLAGS:STRING=-fPIC\ -w
        -DGKLIB_PATH:PATH=${TPSRC}/modmetis-5.1.0/GKlib
        ${TPSRC}/modmetis-5.1.0
)

SET(METIS_LIB metis CACHE FILEPATH "METIS library" FORCE)
MARK_AS_ADVANCED(METIS_LIB)

LINK_DIRECTORIES(${TPDIST}/lib)
INCLUDE_DIRECTORIES(${TPDIST}/include)

IF (WIN32)
    MESSAGE(STATUS "Build Metis: ${TPDIST}/${LIB_DIR}/${METIS_LIB}.dll")
ELSE ()
    MESSAGE(STATUS "Build Metis: ${TPDIST}/${LIB_DIR}/lib${METIS_LIB}.a")
ENDIF()
