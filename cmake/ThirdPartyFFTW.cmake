########################################################################
#
# ThirdParty configuration for Nektar++
#
# FFTW
#
########################################################################

OPTION(NEKTAR_USE_FFTW
    "Use FFTW routines for performing the Fast Fourier Transform." OFF)

IF (NEKTAR_USE_FFTW)
    # Set some common FFTW search paths for the library.
    SET(FFTW_SEARCH_PATHS $ENV{LD_LIBRARY_PATH} $ENV{FFTW_HOME}/lib)
    FIND_LIBRARY(FFTW_LIBRARY NAMES fftw3 fftw3f PATHS ${FFTW_SEARCH_PATHS})

    FIND_PATH(FFTW_INCLUDE_DIR NAMES fftw3.h CACHE FILEPATH 
        "FFTW include directory.")

    IF (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
        SET(BUILD_FFTW OFF)
    ELSE()
        SET(BUILD_FFTW ON)
    ENDIF ()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_FFTW
        "Build FFTW from ThirdParty" ${BUILD_FFTW}
        "NEKTAR_USE_FFTW" OFF)

    IF (THIRDPARTY_BUILD_FFTW)
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            fftw-3.2.2
            URL ${TPURL}/fftw-3.2.2.tar.gz
            URL_MD5 "b616e5c91218cc778b5aa735fefb61ae"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/fftw-3.2.2
            BINARY_DIR ${TPBUILD}/fftw-3.2.2
            TMP_DIR ${TPBUILD}/fftw-3.2.2-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND
                CC=${CMAKE_C_COMPILER}
                ${TPSRC}/fftw-3.2.2/configure
                --prefix=${TPDIST}
                --libdir=${TPDIST}/lib
                --quiet
                --enable-shared
                --disable-dependency-tracking
        )

        SET(FFTW_LIBRARY fftw3 CACHE FILEPATH
            "FFTW library" FORCE)
        SET(FFTW_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "FFTW include" FORCE)

        LINK_DIRECTORIES(${TPDIST}/lib)

        MESSAGE(STATUS "Build FFTW: ${TPDIST}/lib/lib${FFTW_LIBRARY}.so")
        SET(FFTW_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE ()
        ADD_CUSTOM_TARGET(fftw-3.2.2 ALL)
        MESSAGE(STATUS "Found FFTW: ${FFTW_LIBRARY}")
        SET(FFTW_CONFIG_INCLUDE_DIR ${FFTW_INCLUDE_DIR})
    ENDIF()

    # Test if FFTW path is a system path. Only add to include path if not an
    # implicitly defined CXX include path (due to GCC 6.x now providing its own
    # version of some C header files and -isystem reorders include paths).
    GET_FILENAME_COMPONENT(X "${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES}"  ABSOLUTE)
    GET_FILENAME_COMPONENT(Y ${FFTW_INCLUDE_DIR} ABSOLUTE)
    STRING(FIND "${X}" "${Y}" X_FIND)

    IF (X_FIND EQUAL -1)
        INCLUDE_DIRECTORIES(SYSTEM ${FFTW_INCLUDE_DIR})
    ENDIF()

    MARK_AS_ADVANCED(FFTW_LIBRARY)
    MARK_AS_ADVANCED(FFTW_INCLUDE_DIR)
ENDIF( NEKTAR_USE_FFTW )
