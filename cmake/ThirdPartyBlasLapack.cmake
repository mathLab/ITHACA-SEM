########################################################################
#
# ThirdParty configuration for Nektar++
#
# BLAS/LAPACK
#
########################################################################


INCLUDE(FindNativeBlasLapack)

IF(NATIVE_BLAS_LAPACK_FOUND)
    SET(BUILD_BLAS_LAPACK OFF)
ELSE()
    IF(CMAKE_Fortran_COMPILER)
        SET(BUILD_BLAS_LAPACK ON)
    ELSE()
        MESSAGE(SEND_ERROR "No blas installation or fortran compiler found")
    ENDIF()
ENDIF()

OPTION(THIRDPARTY_BUILD_BLAS_LAPACK "Build blas and lapack libraries from ThirdParty."
    ${BUILD_BLAS_LAPACK})

IF(THIRDPARTY_BUILD_BLAS_LAPACK)
    INCLUDE(ExternalProject)

    EXTERNALPROJECT_ADD(
        lapack-3.7.0
        PREFIX ${TPSRC}
        URL http://www.netlib.org/lapack/lapack-3.7.0.tgz
        URL_MD5 "697bb8d67c7d336a0f339cc9dd0fa72f"
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/lapack-3.7.0
        BINARY_DIR ${TPBUILD}/lapack-3.7.0
        TMP_DIR ${TPBUILD}/lapack-3.7.0-tmp
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            -DCMAKE_INSTALL_LIBDIR:PATH=${TPDIST}/lib
            -DBUILD_SHARED_LIBS:STRING=ON
            -DBUILD_TESTING:STRING=OFF
            ${TPSRC}/lapack-3.7.0
        )

    SET(NATIVE_BLAS blas CACHE FILEPATH "BLAS library" FORCE)
    SET(NATIVE_LAPACK lapack CACHE FILEPATH "LAPACK library" FORCE)
    SET(NATIVE_BLAS_LIB_DIR ${TPDIST}/lib CACHE FILEPATH "BLAS library dir" FORCE)
    SET(NATIVE_LAPACK_LIB_DIR ${TPDIST}/lib CACHE FILEPATH "LAPACK library dir" FORCE)
    MARK_AS_ADVANCED(NATIVE_BLAS)
    MARK_AS_ADVANCED(NATIVE_LAPACK)
    MARK_AS_ADVANCED(NATIVE_BLAS_LIB_DIR)
    MARK_AS_ADVANCED(NATIVE_LAPACK_LIB_DIR)

    LINK_DIRECTORIES(${TPDIST}/lib)
    INCLUDE_DIRECTORIES(${TPDIST}/include)

    IF (WIN32)
        MESSAGE(STATUS "Build blas: ${TPDIST}/${LIB_DIR}/${NATIVE_BLAS}.dll")
        MESSAGE(STATUS "Build lapack: ${TPDIST}/${LIB_DIR}/${NATIVE_LAPACK}.dll")
    ELSE ()
        MESSAGE(STATUS "Build blas: ${TPDIST}/${LIB_DIR}/lib${NATIVE_BLAS}.a")
        MESSAGE(STATUS "Build lapack: ${TPDIST}/${LIB_DIR}/${NATIVE_LAPACK}.a")
    ENDIF()
ENDIF()
