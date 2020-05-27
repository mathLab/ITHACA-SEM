########################################################################
#
# ThirdParty configuration for Nektar++
#
# Star CCM i/o
#
########################################################################

OPTION(NEKTAR_USE_CCM "Use CCMIO library for binary Star-CCM+ in NekMesh" OFF)

IF(NEKTAR_USE_CCM)

# First search for system ccmioL installs. Hint /usr/local 
FIND_PATH   (CCMIO_INCLUDE_DIR libccmio/ccmio.h PATHS ${CCM_DIR}/include)
FIND_LIBRARY(CCMIO_LIBRARY NAMES "ccmio" PATHS ${CCM_DIR} PATH_SUFFIXES lib)
FIND_LIBRARY(CCMIO_ADF_LIBRARY NAMES "adf" PATHS ${CCM_DIR} PATH_SUFFIXES lib)

# If we have our library then don't build CCMIO.
IF (CCMIO_INCLUDE_DIR AND CCMIO_LIBRARY AND CCMIO_ADF_LIBRARY)
    SET(BUILD_CCMIO OFF)
ELSE()
    SET(BUILD_CCMIO ON)
ENDIF ()

OPTION(THIRDPARTY_BUILD_CCMIO "Build CCMIO library from ThirdParty." ${BUILD_CCMIO})

IF (THIRDPARTY_BUILD_CCMIO)
    UNSET(PATCH CACHE)
    FIND_PROGRAM(PATCH patch)
    IF(NOT PATCH)
        MESSAGE(FATAL_ERROR
            "'patch' tool for modifying files not found. Cannot build CCM.")
    ENDIF()
    MARK_AS_ADVANCED(PATCH)

    INCLUDE(ExternalProject)
    EXTERNALPROJECT_ADD(
        libccmio-2.6.1
        PREFIX ${TPSRC}
        URL http://visit.ilight.com/svn/visit/trunk/third_party/libccmio-2.6.1.tar.gz
        URL_MD5 f81fbdfb960b1a4f3bcc7feee491efe4
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/libccmio-2.6.1
        BINARY_DIR ${TPBUILD}/libccmio-2.6.1
        TMP_DIR ${TPBUILD}/libccmio-2.6.1-tmp
        PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/cmake/thirdparty-patches/CMakeLists_CCM.txt ${TPSRC}/libccmio-2.6.1/CMakeLists.txt
        COMMAND ${PATCH} -p 0 < ${CMAKE_SOURCE_DIR}/cmake/thirdparty-patches/ccmio-warning.patch
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            ${TPSRC}/libccmio-2.6.1
        )

    THIRDPARTY_LIBRARY(
        CCMIO_LIBRARY STATIC ccmio DESCRIPTION "CCMIO library")
    THIRDPARTY_LIBRARY(
        CCMIO_ADF_LIBRARY STATIC adf DESCRIPTION "CCMIO ADF library")

    INCLUDE_DIRECTORIES(SYSTEM NekMesh ${TPDIST}/include)
    MESSAGE(STATUS "Build CCMIO: ${CCMIO_LIBRARY}")
    SET(CCMIO_CONFIG_INCLUDE_DIR ${TPINC})
ELSE()
    ADD_CUSTOM_TARGET(libccmio-2.6.1 ALL)
    MESSAGE(STATUS "Found CCMIO: ${CCMIO_LIBRARY}")
    SET(CCMIO_CONFIG_INCLUDE_DIR ${CCMIO_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(SYSTEM NekMesh ${CCMIO_INCLUDE_DIR})
    LINK_DIRECTORIES(${CCMIO_LIBRARY_DIR})
ENDIF (THIRDPARTY_BUILD_CCMIO)

SET(CCMIO_LIBRARIES ${CCMIO_LIBRARY} ${CCMIO_ADF_LIBRARY})

MARK_AS_ADVANCED(CCMIO_INCLUDE_DIR)
MARK_AS_ADVANCED(CCMIO_LIBRARY)
MARK_AS_ADVANCED(CCMIO_ADF_LIBRARY)

ENDIF(NEKTAR_USE_CCM)
