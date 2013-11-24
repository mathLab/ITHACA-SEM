OPTION(NEKTAR_USE_MPI "Use MPICH2 for parallelisation." OFF)

CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_GSMPI
    "Build GSMPI if needed" ON
    "NEKTAR_USE_MPI" OFF)

IF( NEKTAR_USE_MPI )
    # First check to see if our compiler has MPI built in to avoid linking libraries etc.
    INCLUDE (CheckIncludeFiles)
    INCLUDE (CheckFunctionExists)
    CHECK_INCLUDE_FILES  (mpi.h    HAVE_MPI_H)
    CHECK_FUNCTION_EXISTS(MPI_Send HAVE_MPI_SEND)

    SET(MPI_BUILTIN OFF CACHE INTERNAL
        "Determines whether MPI is built in")
    IF (NOT "${HAVE_MPI_H}" OR NOT "${HAVE_MPI_SEND}")
        INCLUDE (FindMPI)
        MARK_AS_ADVANCED(MPI_LIBRARY)
        MARK_AS_ADVANCED(MPI_EXTRA_LIBRARY)
        MARK_AS_ADVANCED(file_cmd)
        INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )
        MESSAGE(STATUS "Found MPI: ${MPI_LIBRARY}")
    ELSE()
        SET(MPI_BUILTIN ON)
        MESSAGE(STATUS "Found MPI: built in")
    ENDIF()

    ADD_DEFINITIONS(-DNEKTAR_USE_MPI)

    IF (THIRDPARTY_BUILD_GSMPI)
        EXTERNALPROJECT_ADD(
            gsmpi-1.2
            PREFIX ${TPSRC}
            URL ${TPURL}/gsmpi-1.2.tar.bz2
            URL_MD5 35901be16791bfdeafa9c4d0e06d189b
            DOWNLOAD_DIR ${TPSRC}
            CONFIGURE_COMMAND 
                ${CMAKE_COMMAND}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_BUILD_TYPE:STRING=Debug
                -DCMAKE_INSTALL_PREFIX:PATH=${TPSRC}/dist
                ${TPSRC}/src/gsmpi-1.2
        )
        SET(GSMPI_LIBRARY gsmpi CACHE FILEPATH
            "GSMPI path" FORCE)
        MARK_AS_ADVANCED(GSMPI_LIBRARY)
        SET(XXT_LIBRARY xxt CACHE FILEPATH
            "XXT path" FORCE)
        MARK_AS_ADVANCED(XXT_LIBRARY)
        MESSAGE(STATUS "Build GSMPI: ${TPSRC}/dist/lib/lib${GSMPI_LIBRARY}.a")
        MESSAGE(STATUS "Build XXT: ${TPSRC}/dist/lib/lib${XXT_LIBRARY}.a")
    ELSE (THIRDPARTY_BUILD_GSMPI)
        INCLUDE (FindGSMPI)
        INCLUDE (FindXXT)
    ENDIF (THIRDPARTY_BUILD_GSMPI)
ENDIF( NEKTAR_USE_MPI )

