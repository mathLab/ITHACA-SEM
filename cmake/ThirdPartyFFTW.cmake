# FFTW
SET(NEKTAR_USE_FFTW   OFF CACHE BOOL 
    "Use FFTW routines for performing the Fast Fourier Transform.")
SET(THIRDPARTY_BUILD_FFTW OFF CACHE BOOL
    "Build FFTW from ThirdParty")

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
        INCLUDE_DIRECTORIES(${TPSRC}/dist/include)
    ELSE ()
        INCLUDE (FindFFTW)
        INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
    ENDIF()
    ADD_DEFINITIONS(-DNEKTAR_USING_FFTW)
ENDIF( NEKTAR_USE_FFTW )
