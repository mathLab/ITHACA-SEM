SET(THIRDPARTY_BUILD_ZLIB OFF CACHE BOOL 
    "Build ZLib library")

IF (THIRDPARTY_BUILD_ZLIB)
    # Build the Zlib library separately
    EXTERNALPROJECT_ADD(
        zlib
        PREFIX ${TPSRC}
        URL ${TPURL}/zlib-1.2.3.tar.bz2
        URL_MD5 "dee233bf288ee795ac96a98cc2e369b6"
        DOWNLOAD_DIR ${TPSRC}
        CONFIGURE_COMMAND ./configure --shared --prefix=${TPSRC}/dist
        BUILD_IN_SOURCE 1
    )
    SET(Boost_ZLIB_LIBRARY z)
    SET(Boost_ZLIB_LIBRARY_DEBUG z)
    SET(Boost_ZLIB_LIBRARY_RELEASE z)
ELSE (THIRDPARTY_BUILD_ZLIB)
    # Find a system ZLIB library and use that instead
    FIND_PACKAGE( ZLIB )
    IF (ZLIB_FOUND)
        SET(Boost_ZLIB_LIBRARY ${ZLIB_LIBRARIES})
        SET(Boost_ZLIB_LIBRARY_RELEASE ${ZLIB_LIBRARIES})
        SET(Boost_ZLIB_LIBRARY_DEBUG ${ZLIB_LIBRARIES})
    ENDIF()
ENDIF (THIRDPARTY_BUILD_ZLIB)


