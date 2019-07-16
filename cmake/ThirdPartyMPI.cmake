########################################################################
#
# ThirdParty configuration for Nektar++
#
# MPI
#
########################################################################

OPTION(NEKTAR_USE_MPI "Use MPI for parallelisation." OFF)

CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_GSMPI
    "Build GSMPI if needed" ON
    "NEKTAR_USE_MPI" OFF)

IF( NEKTAR_USE_MPI )
    # First check to see if our compiler has MPI built in to avoid linking
    # libraries etc.
    INCLUDE (CheckIncludeFiles)
    INCLUDE (CheckFunctionExists)
    CHECK_INCLUDE_FILES  (mpi.h    HAVE_MPI_H)
    CHECK_FUNCTION_EXISTS(MPI_Send HAVE_MPI_SEND)

    SET(MPI_BUILTIN OFF CACHE INTERNAL
        "Determines whether MPI is built into the compiler")
    IF (NOT "${HAVE_MPI_H}" OR NOT "${HAVE_MPI_SEND}")
        FIND_PACKAGE(MPI REQUIRED)

        INCLUDE_DIRECTORIES(SYSTEM ${MPI_CXX_INCLUDE_PATH} )
        MESSAGE(STATUS "Found MPI: ${MPI_CXX_LIBRARIES}")
    ELSE()
        SET(MPI_BUILTIN ON)
        MESSAGE(STATUS "Found MPI: built in")
	FIND_PROGRAM(HAVE_APRUN aprun)
	IF (HAVE_APRUN)
	    # Probably on Cray
            SET(MPIEXEC "aprun" CACHE STRING "MPI job launching command")
	    SET(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING
                "MPI job launcher flag to specify number of processes")
	ELSE()
            SET(MPIEXEC "mpirun" CACHE STRING "MPI job launching command")
	    SET(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING
                "MPI job launcher flag to specify number of processes")
	ENDIF()
	MARK_AS_ADVANCED(MPIEXEC)
	MARK_AS_ADVANCED(MPIEXEC_NUMPROC_FLAG)
	UNSET(HAVE_APRUN CACHE)
    ENDIF()

    ADD_DEFINITIONS(-DNEKTAR_USE_MPI)

    IF (THIRDPARTY_BUILD_GSMPI)
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            gsmpi-1.2.1
            URL ${TPURL}/gsmpi-1.2.1_2.tar.bz2
            URL_MD5 3690b4a658324cd5ad17cf8f6115fb5d
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/gsmpi-1.2.1
            BINARY_DIR ${TPBUILD}/gsmpi-1.2.1
            TMP_DIR ${TPBUILD}/gsmpi-1.2.1-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND
                ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_BUILD_TYPE:STRING=Debug
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                ${TPSRC}/gsmpi-1.2.1
        )
        THIRDPARTY_LIBRARY(GSMPI_LIBRARY STATIC gsmpi DESCRIPTION "GSMPI Library")
        THIRDPARTY_LIBRARY(XXT_LIBRARY STATIC xxt DESCRIPTION "XXT Library")
        MARK_AS_ADVANCED(GSMPI_LIBRARY)
        MARK_AS_ADVANCED(XXT_LIBRARY)
        MESSAGE(STATUS "Build GSMPI: ${GSMPI_LIBRARY}")
        MESSAGE(STATUS "Build XXT: ${XXT_LIBRARY}")
    ELSE (THIRDPARTY_BUILD_GSMPI)
        MESSAGE(FATAL_ERROR "Must build GSMPI and XXT")
    ENDIF (THIRDPARTY_BUILD_GSMPI)
ENDIF( NEKTAR_USE_MPI )

