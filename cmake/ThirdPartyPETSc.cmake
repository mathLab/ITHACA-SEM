OPTION(NEKTAR_USE_PETSC
    "Enable PETSc parallel matrix solver support." OFF)

CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_PETSC
    "Build PETSc if needed" OFF
    "NEKTAR_USE_PETSC" OFF)

IF (NEKTAR_USE_PETSC)
    ADD_DEFINITIONS(-DNEKTAR_USING_PETSC)

    IF (THIRDPARTY_BUILD_PETSC)
        INCLUDE(ExternalProject)

        SET(PETSC_C_COMPILER "${CMAKE_C_COMPILER}")
        SET(PETSC_CXX_COMPILER "${CMAKE_CXX_COMPILER}")

        IF (NEKTAR_USE_MPI)
            IF (NOT MPI_BUILTIN)
                SET(PETSC_C_COMPILER "${MPI_C_COMPILER}")
                SET(PETSC_CXX_COMPILER "${MPI_CXX_COMPILER}")
            ENDIF (NOT MPI_BUILTIN)
        ELSE (NEKTAR_USE_MPI)
            SET(PETSC_NO_MPI "--with-mpi=0")
        ENDIF (NEKTAR_USE_MPI)

        EXTERNALPROJECT_ADD(
            petsc-3.5.1
            PREFIX ${TPSRC}
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/petsc-3.5.1
            TMP_DIR ${TPBUILD}/petsc-3.5.1-tmp
            INSTALL_DIR ${TPDIST}
            BINARY_DIR ${TPBUILD}/petsc-3.5.1
            URL http://www.nektar.info/thirdparty/petsc-lite-3.5.1.tar.gz
            URL_MD5 "539b3bdb627407b7e4e9e830fd5ccf43"
            CONFIGURE_COMMAND ./configure
                --with-cc=${PETSC_C_COMPILER}
                --with-cxx=${PETSC_CXX_COMPILER}
                --with-shared-libraries=0
                --with-pic=1
                --with-x=0
                --prefix=${TPDIST}
                --with-petsc-arch=c-opt
                --with-fc=0
                ${PETSC_NO_MPI}
                )

        INCLUDE_DIRECTORIES(${TPDIST}/include)
        SET(PETSC_LIBRARIES "${TPDIST}/lib/libpetsc.a")
        MESSAGE(STATUS "Build PETSc: ${PETSC_LIBRARIES}")
    ELSE (THIRDPARTY_BUILD_PETSC)
        INCLUDE(FindPETSc)
        IF (NOT PETSC_FOUND)
            MESSAGE(FATAL_ERROR "Could not find PETSc")
        ELSE (NOT PETSC_FOUND)
            MESSAGE(STATUS "Found PETSc: ${PETSC_LIBRARIES}")
        ENDIF (NOT PETSC_FOUND)
        INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
    ENDIF (THIRDPARTY_BUILD_PETSC)
ENDIF( NEKTAR_USE_PETSC )

INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
