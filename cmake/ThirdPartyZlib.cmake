########################################################################
#
# ThirdParty configuration for Nektar++
#
# ZLib
#
########################################################################

# Attempt to identify Macports libraries, if they exist. This prevents
# cmake warnings later on.
SET(ZLIB_ROOT /opt/local/)


# Find a system ZLIB library
# If not found enable the THIRDPARTY_BUILD_ZLIB option
FIND_PACKAGE( ZLIB )
IF (ZLIB_FOUND AND NOT ZLIB_VERSION_PATCH LESS 7 )
    SET(ZLIB_LIBRARY ${ZLIB_LIBRARIES})
    SET(ZLIB_LIBRARY_DEBUG ${ZLIB_LIBRARIES})
    SET(ZLIB_LIBRARY_RELEASE ${ZLIB_LIBRARIES})
    MESSAGE(STATUS "Found Zlib library: ${ZLIB_LIBRARY}")
    OPTION(THIRDPARTY_BUILD_ZLIB "Build ZLib library" OFF)
ELSE ()
    OPTION(THIRDPARTY_BUILD_ZLIB "Build ZLib library" ON)
ENDIF()


# If we or the user
IF (THIRDPARTY_BUILD_ZLIB)
    MESSAGE(STATUS "Will build Zlib 1.2.7")
    # Build the Zlib library separately
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
    IF (WIN32)
        SET(ZLIB_LIBRARY zlib)
        SET(ZLIB_LIBRARY_DEBUG zlibd)
        SET(ZLIB_LIBRARY_RELEASE zlib)
    ELSE (WIN32)
        SET(ZLIB_LIBRARY z)
        SET(ZLIB_LIBRARY_DEBUG z)
        SET(ZLIB_LIBRARY_RELEASE z)
    ENDIF (WIN32)
ELSE (THIRDPARTY_BUILD_ZLIB)
    ADD_CUSTOM_TARGET(zlib-1.2.7 ALL)
ENDIF (THIRDPARTY_BUILD_ZLIB)

