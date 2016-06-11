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
    # First check to see if our compiler has MPI built in to avoid linking libraries etc.
    INCLUDE (CheckIncludeFiles)
    INCLUDE (CheckFunctionExists)
    CHECK_INCLUDE_FILES  (mpi.h    HAVE_MPI_H)
    CHECK_FUNCTION_EXISTS(MPI_Send HAVE_MPI_SEND)

    SET(BUILD_MPI OFF)
    SET(MPI_BUILTIN OFF CACHE INTERNAL
        "Determines whether MPI is built into the compiler")
    IF (NOT "${HAVE_MPI_H}" OR NOT "${HAVE_MPI_SEND}")
        INCLUDE (FindMPI)
        MARK_AS_ADVANCED(MPI_LIBRARY)
        MARK_AS_ADVANCED(MPI_EXTRA_LIBRARY)
        MARK_AS_ADVANCED(file_cmd)
        INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )

        IF (NOT MPI_CXX_FOUND OR NOT MPI_C_FOUND)
            SET(BUILD_MPI ON)
        ELSE()
            MESSAGE(STATUS "Found MPI: ${MPI_LIBRARY}")
        ENDIF()
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

    IF (BUILD_MPI)
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            mpich-3.2
            PREFIX ${TPSRC}
            URL http://ae-nektar.ae.ic.ac.uk/~dmoxey/mpich-3.2.tar.gz
            URL_MD5 f414cfa77099cd1fa1a5ae4e22db508a
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/mpich-3.2
            BINARY_DIR ${TPBUILD}/mpich-3.2
            TMP_DIR ${TPBUILD}/mpich-3.2
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ./configure --prefix=${TPDIST} --disable-f77 --disable-fc --disable-mpe
            )
    ENDIF()

    IF (THIRDPARTY_BUILD_GSMPI)
        EXTERNALPROJECT_ADD(
            gsmpi-1.2
            URL ${TPURL}/gsmpi-1.2.tar.bz2
            URL_MD5 35901be16791bfdeafa9c4d0e06d189b
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/gsmpi-1.2
            BINARY_DIR ${TPBUILD}/gsmpi-1.2
            TMP_DIR ${TPBUILD}/gsmpi-1.2-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND 
                ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_BUILD_TYPE:STRING=Debug
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                ${TPSRC}/gsmpi-1.2
        )
        SET(GSMPI_LIBRARY gsmpi CACHE FILEPATH
            "GSMPI path" FORCE)
        MARK_AS_ADVANCED(GSMPI_LIBRARY)
        SET(XXT_LIBRARY xxt CACHE FILEPATH
            "XXT path" FORCE)
        MARK_AS_ADVANCED(XXT_LIBRARY)
        THIRDPARTY_STATIC_LIBNAME(GSMPI_LIBRARY)
        THIRDPARTY_STATIC_LIBNAME(XXT_LIBRARY)
        MESSAGE(STATUS "Build GSMPI: ${GSMPI_LIBRARY}")
        MESSAGE(STATUS "Build XXT: ${XXT_LIBRARY}")
    ELSE (THIRDPARTY_BUILD_GSMPI)
        ADD_CUSTOM_TARGET(gsmpi-1.2 ALL)
        INCLUDE (FindGSMPI)
        INCLUDE (FindXXT)
    ENDIF (THIRDPARTY_BUILD_GSMPI)
ENDIF( NEKTAR_USE_MPI )

