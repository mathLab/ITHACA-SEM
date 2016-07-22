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

    SET(BUILD_MPI OFF)
    SET(MPI_BUILTIN OFF CACHE INTERNAL
        "Determines whether MPI is built into the compiler")

    IF (NOT "${HAVE_MPI_H}" OR NOT "${HAVE_MPI_SEND}")
        # No in-built MPI: try to find it on the system instead.
        IF (NOT THIRDPARTY_BUILD_MPI)
            INCLUDE (FindMPI)

            IF (NOT MPI_CXX_FOUND)
                # No MPI at all: we have to build it
                SET(BUILD_MPI ON)
            ELSE()
                MARK_AS_ADVANCED(file_cmd)
                INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
                MESSAGE(STATUS "Found MPI: ${MPI_CXX_LIBRARIES}")
            ENDIF()
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

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_MPI
        "Build MPI library from ThirdParty" ${BUILD_MPI}
        "NEKTAR_USE_MPI" OFF)

    IF (THIRDPARTY_BUILD_MPI)
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            openmpi-1.10.3
            PREFIX ${TPSRC}
            URL http://ae-nektar.ae.ic.ac.uk/~dmoxey/openmpi-1.10.3.tar.gz
            URL_MD5 7d384d6eb454d3a621f932058a822b14
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/openmpi-1.10.3
            BINARY_DIR ${TPBUILD}/openmpi-1.10.3
            TMP_DIR ${TPBUILD}/openmpi-1.10.3
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${TPSRC}/openmpi-1.10.3/configure --prefix=${TPDIST} --disable-mpi-fortran --disable-static --disable-vt
            )

        SET_TARGET_PROPERTIES(openmpi-1.10.3 PROPERTIES EXCLUDE_FROM_ALL 1)
        THIRDPARTY_LIBRARY(MPI_CXX_LIBRARIES SHARED mpi mpi_cxx
            DESCRIPTION "MPI C++ libraries")
        SET(MPI_CXX_COMPILE_FLAGS "" CACHE STRING "MPI compiler flags" FORCE)
        SET(MPI_CXX_INCLUDE_PATH ${TPDIST}/include CACHE FILEPATH "MPI include path" FORCE)
        SET(MPI_CXX_COMPILER ${TPDIST}/bin/mpicxx CACHE FILEPATH "MPI C++ compiler" FORCE)
        SET(MPI_C_COMPILER ${TPDIST}/bin/mpicc CACHE FILEPATH "MPI C compiler" FORCE)

        SET(OMPI_EXES mpiexec mpirun orted orterun)

        FOREACH(exe ${OMPI_EXES})
            LIST(APPEND OMPI_EXES2 "${TPDIST}/bin/${exe}")
        ENDFOREACH()

        INSTALL(PROGRAMS ${OMPI_EXES2}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/${NEKTAR_BIN_DIR}
            COMPONENT ThirdParty)

        MESSAGE(STATUS "Build MPI: ${MPI_CXX_LIBRARIES}")
    ELSE ()
        ADD_CUSTOM_TARGET(openmpi-1.10.3 ALL)
        SET(MPI_CONFIG_INCLUDE_DIR ${MPI_CXX_INCLUDE_PATH})
    ENDIF()

    EXTERNALPROJECT_ADD(
        gsmpi-1.2.1
        URL ${TPURL}/gsmpi-1.2.1.tar.bz2
        URL_MD5 18dcb4cd1dcc7876173465c404b1142d
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


    THIRDPARTY_LIBRARY(GSMPI_LIBRARY STATIC gsmpi
        DESCRIPTION "GSMPI library")
    THIRDPARTY_LIBRARY(XXT_LIBRARY STATIC xxt
        DESCRIPTION "XXT library")
    MARK_AS_ADVANCED(GSMPI_LIBRARY)
    MARK_AS_ADVANCED(XXT_LIBRARY)
    MESSAGE(STATUS "Build GSMPI: ${GSMPI_LIBRARY}")
    MESSAGE(STATUS "Build XXT: ${XXT_LIBRARY}")

    ADD_DEPENDENCIES(gsmpi-1.2.1 openmpi-1.10.3)
ENDIF(NEKTAR_USE_MPI)
