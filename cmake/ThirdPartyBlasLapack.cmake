########################################################################
#
# ThirdParty configuration for Nektar++
#
# BLAS/LAPACK
#
########################################################################

# Test for if BLAS is in compiler.
INCLUDE(CheckFunctionExists)

CHECK_FUNCTION_EXISTS(dgemm_ HAVE_DGEMM)

IF (HAVE_DGEMM)
    MESSAGE(STATUS "Found BLAS/LAPACK: built in")
    SET(BLAS_LAPACK_BUILTIN ON)
ELSE()
    SET(BLAS_LAPACK_BUILTIN OFF)

    # BLAS Support
    CMAKE_DEPENDENT_OPTION(NEKTAR_USE_SYSTEM_BLAS_LAPACK
        "Use the system provided blas and lapack libraries" ON
        "UNIX; NOT APPLE; NOT NEKTAR_USE_OPENBLAS; NOT NEKTAR_USE_MKL; NOT NEKTAR_USE_ACML; NOT NEKTAR_USE_ACCELERATE_FRAMEWORK" OFF)
    CMAKE_DEPENDENT_OPTION(NEKTAR_USE_OPENBLAS
        "Use OpenBLAS library as a substitute to native BLAS." OFF
        "NOT NEKTAR_USE_SYSTEM_BLAS_LAPACK" OFF)
    CMAKE_DEPENDENT_OPTION(NEKTAR_USE_ACML
        "Use the AMD Core Math Library (ACML) for BLAS and Lapack support." OFF
        "NOT NEKTAR_USE_SYSTEM_BLAS_LAPACK" OFF)
    CMAKE_DEPENDENT_OPTION(NEKTAR_USE_MKL
        "Use the Intel Math Kernel Library (MKL) for BLAS and Lapack support." OFF
        "NOT NEKTAR_USE_SYSTEM_BLAS_LAPACK" OFF)
    CMAKE_DEPENDENT_OPTION(NEKTAR_USE_ACCELERATE_FRAMEWORK
        "Use the Mac Accelerate Framework for BLAS and Lapack support." ON
        "NOT NEKTAR_USE_SYSTEM_BLAS_LAPACK; APPLE" OFF)
    CMAKE_DEPENDENT_OPTION(NEKTAR_USE_WIN32_LAPACK
        "Use Win32 Lapack provided with the Third Party Distribution."
        ON "NOT NEKTAR_USE_SYSTEM_BLAS_LAPACK; WIN32" OFF)

    INCLUDE(FindBlasLapack)

    IF(LAPACK_FOUND)
        SET(BUILD_BLAS_LAPACK OFF)
    ELSE()
        IF(CMAKE_Fortran_COMPILER AND NATIVE_BLAS_LAPACK_FIND_REQUIRED)
            SET(BUILD_BLAS_LAPACK ON)
        ELSE()
            MESSAGE(SEND_ERROR "No blas/lapack installation found and cannot self build")
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

        SET(BLAS_LAPACK blas lapack)

        LINK_DIRECTORIES(${TPDIST}/lib)
        INCLUDE_DIRECTORIES(${TPDIST}/include)

        IF (WIN32)
            MESSAGE(STATUS "Build blas: ${TPDIST}/${LIB_DIR}/libblas.dll")
            MESSAGE(STATUS "Build lapack: ${TPDIST}/${LIB_DIR}/liblapack.dll")
        ELSE ()
            MESSAGE(STATUS "Build blas: ${TPDIST}/${LIB_DIR}/libblas.a")
            MESSAGE(STATUS "Build lapack: ${TPDIST}/${LIB_DIR}/liblapack.a")
        ENDIF()
    ENDIF()
ENDIF()