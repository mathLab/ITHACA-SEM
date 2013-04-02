OPTION(THIRDPARTY_BUILD_ZLIB "Build ZLib library" OFF)

IF (THIRDPARTY_BUILD_ZLIB)
    # Build the Zlib library separately
    EXTERNALPROJECT_ADD(
        zlib-1.2.7
        PREFIX ${TPSRC}
        URL ${TPURL}/zlib-1.2.7.tar.gz
        URL_MD5 "60df6a37c56e7c1366cca812414f7b85"
        DOWNLOAD_DIR ${TPSRC}
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX:PATH=${TPSRC}/dist -DCMAKE_C_FLAGS:STRING=-fPIC ${TPSRC}/src/zlib-1.2.7
    )
    SET(ZLIB_LIBRARY z)
    SET(ZLIB_LIBRARY_DEBUG z)
    SET(ZLIB_LIBRARY_RELEASE z)
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
  
  #  IF (ZLIB_LIBRARY)
  #  SET(Boost_ZLIB_LIBRARY ${ZLIB_LIBRARY})
  #  SET(Boost_ZLIB_LIBRARY_RELEASE ${ZLIB_LIBRARY})
  #  SET(Boost_ZLIB_LIBRARY_DEBUG ${ZLIB_LIBRARY})
  #ENDIF()
  
  IF (ZLIB_FOUND)
    MESSAGE(STATUS "Found Zlib library: ${ZLIB_LIBRARY}")
  ENDIF ()
  
ENDIF (THIRDPARTY_BUILD_ZLIB)


