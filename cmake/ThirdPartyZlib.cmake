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
        PREFIX ${TPSRC}
        URL ${TPURL}/zlib-1.2.7.tar.gz
        URL_MD5 "4a162e0f643232e7e278d59a0603ceb0"
        DOWNLOAD_DIR ${TPSRC}
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX:PATH=${TPSRC}/dist -DCMAKE_C_FLAGS:STRING=-fPIC ${TPSRC}/src/zlib-1.2.7
    )
    IF (WIN32)
        SET(ZLIB_LIBRARY zlib)
        SET(ZLIB_LIBRARY_DEBUG zlib)
        SET(ZLIB_LIBRARY_RELEASE zlib)
    ELSE (WIN32)
        SET(ZLIB_LIBRARY z)
        SET(ZLIB_LIBRARY_DEBUG z)
        SET(ZLIB_LIBRARY_RELEASE z)
    ENDIF (WIN32)
ENDIF (THIRDPARTY_BUILD_ZLIB)

