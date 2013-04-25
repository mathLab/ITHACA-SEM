# FFTW
OPTION(NEKTAR_USE_FFTW
    "Use FFTW routines for performing the Fast Fourier Transform." OFF)

CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_FFTW
    "Build FFTW from ThirdParty" OFF
    "NEKTAR_USE_FFTW" OFF)

IF( NEKTAR_USE_FFTW )
    IF (THIRDPARTY_BUILD_FFTW)
        INCLUDE(ExternalProject)
        
        EXTERNALPROJECT_ADD(
            fftw-3.2.2
            PREFIX ${TPSRC}
            URL ${TPURL}/fftw-3.2.2.tar.gz
            URL_MD5 "b616e5c91218cc778b5aa735fefb61ae"
            DOWNLOAD_DIR ${TPSRC}
            CONFIGURE_COMMAND ${TPSRC}/src/fftw-3.2.2/configure --prefix=${TPSRC}/dist --quiet --enable-shared --disable-dependency-tracking
        )
        SET(FFTW_LIB fftw3)
        MARK_AS_ADVANCED(FFTW_LIB)
        INCLUDE_DIRECTORIES(${TPSRC}/dist/include)
        MESSAGE(STATUS "Build FFTW: ${TPSRC}/dist/lib/lib${FFTW_LIB}.so")
    ELSE ()
        INCLUDE (FindFFTW)
        INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
    ENDIF()
ENDIF( NEKTAR_USE_FFTW )
