########################################################################
#
# ThirdParty configuration for Nektar++
#
# ZLib
#
########################################################################

# Attempt to identify Macports libraries, if they exist and we didn't override
# ZLIB_ROOT on the command line or in ccmake. This prevents cmake warnings later
# on.
IF (NOT DEFINED ZLIB_ROOT)
    SET(ZLIB_ROOT /opt/local/)
ENDIF()

# Find a system ZLIB library. If not found enable the THIRDPARTY_BUILD_ZLIB
# option.
SET(ZLIB_FIND_QUIETLY ON)
FIND_PACKAGE(ZLIB)
IF (ZLIB_FOUND AND NOT ZLIB_VERSION_PATCH LESS 7)
    SET(ZLIB_LIBRARY ${ZLIB_LIBRARIES} CACHE FILEPATH
        "Zlib library" FORCE)
    SET(ZLIB_LIBRARY_DEBUG ${ZLIB_LIBRARIES} CACHE FILEPATH
        "Zlib library" FORCE)
    SET(ZLIB_LIBRARY_RELEASE ${ZLIB_LIBRARIES} CACHE FILEPATH
        "Zlib library" FORCE)
    MARK_AS_ADVANCED(ZLIB_LIBRARY ZLIB_LIBRARY_DEBUG ZLIB_LIBRARY_RELEASE)
    SET(ZLIB_INCLUDE_DIR ${ZLIB_INCLUDE_DIRS} CACHE PATH
        "Zlib include" FORCE)
    SET(BUILD_ZLIB OFF)
ELSE ()
    SET(BUILD_ZLIB ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_ZLIB "Build ZLib library" ${BUILD_ZLIB})

# If we or the user
IF (THIRDPARTY_BUILD_ZLIB)
    INCLUDE(ExternalProject)
    EXTERNALPROJECT_ADD(
        zlib-1.2.7
        URL ${TPURL}/zlib-1.2.7.tar.gz
        URL_MD5 "4a162e0f643232e7e278d59a0603ceb0"
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/zlib-1.2.7
        BINARY_DIR ${TPBUILD}/zlib-1.2.7
        TMP_DIR ${TPBUILD}/zlib-1.2.7-tmp
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            -DCMAKE_C_FLAGS:STRING=-fPIC
            ${TPSRC}/zlib-1.2.7
        )

    IF (APPLE)
        EXTERNALPROJECT_ADD_STEP(zlib-1.2.7 patch-install-path
            COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id ${CMAKE_INSTALL_PREFIX}/${NEKTAR_LIB_DIR}/libz.1.2.7.dylib ${TPDIST}/lib/libz.1.2.7.dylib
            DEPENDEES install)
    ENDIF ()

    IF (WIN32)
        THIRDPARTY_LIBRARY(ZLIB_LIBRARY SHARED zlib
            DESCRIPTION "Zlib library" BINDIR)
        THIRDPARTY_LIBRARY(ZLIB_LIBRARY_DEBUG SHARED zlibd
            DESCRIPTION "Zlib library" BINDIR)
        THIRDPARTY_LIBRARY(ZLIB_LIBRARY_RELEASE SHARED zlib
            DESCRIPTION "Zlib library" BINDIR)
    ELSE ()
        THIRDPARTY_LIBRARY(ZLIB_LIBRARY SHARED z
            DESCRIPTION "Zlib library")
        THIRDPARTY_LIBRARY(ZLIB_LIBRARY_DEBUG SHARED z
            DESCRIPTION "Zlib library")
        THIRDPARTY_LIBRARY(ZLIB_LIBRARY_RELEASE SHARED z
            DESCRIPTION "Zlib library")
    ENDIF ()

    MESSAGE(STATUS "Build Zlib: ${ZLIB_LIBRARY}")
    SET(ZLIB_INCLUDE_DIR ${TPDIST}/include CACHE PATH "Zlib include" FORCE)
    SET(ZLIB_CONFIG_INCLUDE_DIR ${TPINC})
ELSE (THIRDPARTY_BUILD_ZLIB)
    ADD_CUSTOM_TARGET(zlib-1.2.7 ALL)
    MESSAGE(STATUS "Found Zlib: ${ZLIB_LIBRARY} (version ${ZLIB_VERSION_STRING})")
    SET(ZLIB_CONFIG_INCLUDE_DIR ${ZLIB_INCLUDE_DIR})
ENDIF (THIRDPARTY_BUILD_ZLIB)

MARK_AS_ADVANCED(ZLIB_LIBRARY)
MARK_AS_ADVANCED(ZLIB_LIBRARY_DEBUG)
MARK_AS_ADVANCED(ZLIB_LIBRARY_RELEASE)
MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
