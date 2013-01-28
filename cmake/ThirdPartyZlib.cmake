OPTION(THIRDPARTY_BUILD_ZLIB "Build ZLib library" OFF)

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
  # Attempt to identify Macports libraries, if they exist. This prevents
  # cmake warnings later on.
  IF (APPLE)
    FIND_LIBRARY(ZLIB_LIBRARY z PATHS /opt/local/lib NO_DEFAULT_PATH)
  ENDIF()
  
  # Find a system ZLIB library and use that instead
  IF (ZLIB_LIBRARY-NOTFOUND OR NOT APPLE) 
    FIND_PACKAGE( ZLIB )
    IF (ZLIB_FOUND)
      SET(ZLIB_LIBRARY ${ZLIB_LIBRARIES})
    ENDIF()
  ENDIF()
  
  IF (ZLIB_LIBRARY)
    SET(Boost_ZLIB_LIBRARY ${ZLIB_LIBRARY})
    SET(Boost_ZLIB_LIBRARY_RELEASE ${ZLIB_LIBRARY})
    SET(Boost_ZLIB_LIBRARY_DEBUG ${ZLIB_LIBRARY})
  ENDIF()
ENDIF (THIRDPARTY_BUILD_ZLIB)


