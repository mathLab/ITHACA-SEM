########################################################################
#
# ThirdParty configuration for Nektar++
#
# ZLib
#
########################################################################

# Find a system ZLIB library. If not found enable the THIRDPARTY_BUILD_ZLIB
# option.
FIND_PACKAGE(ZLIB QUIET)
IF (ZLIB_FOUND AND NOT ZLIB_VERSION_PATCH LESS 7)
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
        THIRDPARTY_LIBRARY(ZLIB_LIBRARIES STATIC zlib DESCRIPTION "Zlib library")
        THIRDPARTY_LIBRARY(ZLIB_LIBRARIES_DEBUG STATIC zlibd DESCRIPTION "Zlib library")
    ELSE ()
        THIRDPARTY_LIBRARY(ZLIB_LIBRARIES SHARED z DESCRIPTION "Zlib library")
        THIRDPARTY_LIBRARY(ZLIB_LIBRARIES_DEBUG SHARED z DESCRIPTION "Zlib library")
    ENDIF ()

    MESSAGE(STATUS "Build Zlib: ")
    MESSAGE(STATUS " -- Optimized: ${ZLIB_LIBRARIES}")
    MESSAGE(STATUS " -- Debug:     ${ZLIB_LIBRARIES_DEBUG}")
    SET(ZLIB_INCLUDE_DIR ${TPDIST}/include CACHE PATH "Zlib include" FORCE)
    SET(ZLIB_CONFIG_INCLUDE_DIR ${TPINC})
ELSE (THIRDPARTY_BUILD_ZLIB)
    ADD_CUSTOM_TARGET(zlib-1.2.7 ALL)
    MESSAGE(STATUS "Found Zlib: ${ZLIB_LIBRARIES} (version ${ZLIB_VERSION_STRING})")

    # We use the found library also for debug builds.
    SET(ZLIB_LIBRARIES_DEBUG ${ZLIB_LIBRARIES})

    SET(ZLIB_CONFIG_INCLUDE_DIR ${ZLIB_INCLUDE_DIRS})
ENDIF (THIRDPARTY_BUILD_ZLIB)

INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
