########################################################################
#
# ThirdParty configuration for Nektar++
#
# PETSc
#
########################################################################

OPTION(NEKTAR_USE_PETSC
    "Enable PETSc parallel matrix solver support." OFF)

IF (NEKTAR_USE_PETSC)
    SET(PETSC_FIND_QUIETLY ON)
    INCLUDE(FindPETSc)

    IF (PETSC_FOUND)
        SET(BUILD_PETSC OFF)
    ELSE()
        SET(BUILD_PETSC ON)
    ENDIF()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_PETSC
        "Build PETSc if needed" ${BUILD_PETSC}
        "NEKTAR_USE_PETSC" OFF)

    IF (THIRDPARTY_BUILD_PETSC)
        INCLUDE(ExternalProject)

        FIND_PACKAGE(PythonInterp 2 REQUIRED)

        SET(PETSC_C_COMPILER   "${CMAKE_C_COMPILER}")
        SET(PETSC_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
        SET(PETSC_Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")

        IF (NEKTAR_USE_MPI)
            IF (NOT MPI_BUILTIN)
                SET(PETSC_C_COMPILER   "${MPI_C_COMPILER}")
                SET(PETSC_CXX_COMPILER "${MPI_CXX_COMPILER}")
                SET(PETSC_Fortran_COMPILER "${MPI_Fortran_COMPILER}")
            ENDIF (NOT MPI_BUILTIN)
        ELSE (NEKTAR_USE_MPI)
            SET(PETSC_NO_MPI "--with-mpi=0")
        ENDIF (NEKTAR_USE_MPI)

        IF(CMAKE_Fortran_COMPILER)
            IF(NEKTAR_USE_MPI AND NOT MPI_Fortran_COMPILER)
                MESSAGE(ERROR "MPI_Fortran_COMPILER not set")
            ENDIF()
            # we use a MUMPS build in ordering here, in the future it might make
            # sense to hook it up with metis/scotch since this MIGHT be faster
            SET(PETSC_MUMPS --download-scalapack --download-mumps)
            SET(PETSC_DEPS "")

            IF( NEKTAR_USE_BLAS_LAPACK )
                IF( NEKTAR_USE_MKL AND MKL_FOUND )
                    SET(PETSC_MUMPS ${PETSC_MUMPS} --with-blas-lapack-dir=${MKL_LIB_DIR})
                ELSEIF( NEKTAR_USE_WIN32_LAPACK )
                    SET(PETSC_MUMPS ${PETSC_MUMPS} --with-blas-lapack-dir=${LAPACK_DIR})
                ELSEIF( NEKTAR_USE_SYSTEM_BLAS_LAPACK )
                    SET(PETSC_MUMPS ${PETSC_MUMPS} --with-blas-lapack-dir=${NATIVE_LAPACK_LIB_DIR})
                    IF(THIRDPARTY_BUILD_BLAS_LAPACK)
                        SET(PETSC_DEPS ${PETSC_DEPS} lapack-3.7.0)
                    ENDIF()
                ELSE()
                    MESSAGE(STATUS "No suitable blas/lapack found, downloading")
                    SET(PETSC_MUMPS ${PETSC_MUMPS} --download-fblaslapack)
                ENDIF()
            ELSE()
                MESSAGE(STATUS "No suitable blas/lapack found, downloading")
                SET(PETSC_MUMPS ${PETSC_MUMPS} --download-fblaslapack)
            ENDIF()

        ELSE()
            MESSAGE(WARNING "No Fortran support. Building PETSc without MUMPS support")
            SET(PETSC_Fortran_COMPILER "0")
        ENDIF()

        EXTERNALPROJECT_ADD(
            petsc-3.7.2
            DEPENDS ${PETSC_DEPS}
            PREFIX ${TPSRC}
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/petsc-3.7.2
            TMP_DIR ${TPBUILD}/petsc-3.7.2-tmp
            INSTALL_DIR ${TPDIST}
            BINARY_DIR ${TPBUILD}/petsc-3.7.2
            URL http://www.nektar.info/thirdparty/petsc-lite-3.7.2.tar.gz
            URL_MD5 "26c2ff8eaaa9e49aea063f839f5daa7e"
            CONFIGURE_COMMAND
                OMPI_FC=${CMAKE_Fortran_COMPILER}
                OMPI_CC=${CMAKE_C_COMPILER}
                OMPI_CXX=${CMAKE_CXX_COMPILER}
                ${PYTHON_EXECUTABLE} ./configure
                ./configure
                --with-fc=${PETSC_Fortran_COMPILER}
                --with-cc=${PETSC_C_COMPILER}
                --with-cxx=${PETSC_CXX_COMPILER}
                --with-shared-libraries=1
                --with-pic=1
                --with-x=0
                --with-ssl=0
                --prefix=${TPDIST}
                --with-petsc-arch=c-opt
                ${PETSC_MUMPS}
                ${PETSC_NO_MPI}
            BUILD_COMMAND MAKEFLAGS= make)

        SET(PETSC_LIBRARIES petsc CACHE FILEPATH
            "PETSc library" FORCE)
        SET(PETSC_INCLUDES ${TPDIST}/include CACHE FILEPATH
            "PETSc includes" FORCE)

        LINK_DIRECTORIES(${TPDIST}/lib)
        MESSAGE(STATUS "Build PETSc: ${TPDIST}/${LIB_DIR}/lib${PETSC_LIBRARIES}.so")
        SET(PETSC_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE (THIRDPARTY_BUILD_PETSC)
        INCLUDE(FindPETSc)
        IF (NOT PETSC_FOUND)
            MESSAGE(FATAL_ERROR "Could not find PETSc")
        ELSE (NOT PETSC_FOUND)
            MESSAGE(STATUS "Found PETSc: ${PETSC_LIBRARIES}")
        ENDIF (NOT PETSC_FOUND)
        SET(PETSC_CONFIG_INCLUDE_DIR ${PETSC_INCLUDES})
        INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
        ADD_CUSTOM_TARGET(petsc-3.7.2 ALL)
    ENDIF (THIRDPARTY_BUILD_PETSC)

    ADD_DEFINITIONS(-DNEKTAR_USING_PETSC)
    INCLUDE_DIRECTORIES(SYSTEM ${PETSC_INCLUDES})
    IF (NOT NEKTAR_USE_MPI)
        INCLUDE_DIRECTORIES(SYSTEM ${PETSC_INCLUDES}/petsc/mpiuni)
    ENDIF (NOT NEKTAR_USE_MPI)

    MARK_AS_ADVANCED(PETSC_CURRENT PETSC_DIR PETSC_LIBRARIES PETSC_INCLUDES)
ENDIF( NEKTAR_USE_PETSC )
