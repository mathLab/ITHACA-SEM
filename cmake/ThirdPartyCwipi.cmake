########################################################################
#
# ThirdParty configuration for Nektar++
#
# CWIPI
#
########################################################################

OPTION(NEKTAR_USE_CWIPI
    "Use CWIPI for Coupling." OFF)

IF ( NEKTAR_USE_CWIPI )
    set(CWIPI_DIR "/usr/local/cwipi" CACHE PATH "CWIPI base path")
    find_package(Cwipi)
    include_directories(${CWIPI_INCLUDE_DIRS})
    ADD_DEFINITIONS(-DNEKTAR_USE_CWIPI)
ENDIF ( NEKTAR_USE_CWIPI )

IF (NEKTAR_USE_CWIPI)
    # Set some common CWIPI search paths.
    SET(CWIPI_SEARCH_PATHS $ENV{LD_LIBRARY_PATH} $ENV{CWIPI_HOME}/lib)
    FIND_LIBRARY(CWIPI_LIBRARY NAMES cwipi fvmc bftc PATHS ${CWIPI_SEARCH_PATHS})

    IF (CWIPI_LIBRARY)
        GET_FILENAME_COMPONENT(CWIPI_PATH ${CWIPI_LIBRARY} PATH)
        SET(CWIPI_INCLUDE_DIR ${CWIPI_PATH}/../include CACHE FILEPATH "CWIPI include directory.")
        SET(BUILD_CWIPI OFF)
    ELSE()
        SET(BUILD_CWIPI ON)
    ENDIF ()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_CWIPI
        "Build CWIPI from ThirdParty" ${BUILD_CWIPI}
        "NEKTAR_USE_CWIPI" OFF)

    IF (THIRDPARTY_BUILD_CWIPI)
        INCLUDE(ExternalProject)

        get_filename_component(MPI_LIBRARY_PATH ${MPI_LIBRARY} DIRECTORY)
        get_filename_component(MPI_BIN_PATH ${MPI_C_COMPILER} DIRECTORY)
        get_filename_component(MPI_INC_PATH ${MPI_LIBRARY_PATH} DIRECTORY)
        set(MPI_INC_PATH "${MPI_INC_PATH}/include")

        EXTERNALPROJECT_ADD(
            cwipi-0.8.2
            URL "http://sites.onera.fr/cwipi/sites/sites.onera.fr.cwipi/files/u4/cwipi-0.8.2.tgz"
            URL_MD5 "cd28dbea20a08d71f5ff4b4770867268"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/cwipi-0.8.2
            BINARY_DIR ${TPBUILD}/cwipi-0.8.2
            TMP_DIR ${TPBUILD}/cwipi-0.8.2-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${TPSRC}/cwipi-0.8.2/configure
                CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
                --with-mpi-include=${MPI_INC_PATH}
                --with-mpi-lib=${MPI_LIBRARY_PATH}
                --with-mpi-bin=${MPI_BIN_PATH}
                --prefix=${TPDIST}
                --libdir=${TPDIST}/lib
                --quiet
        )

        SET(CWIPI_LIBRARY cwipi CACHE FILEPATH
            "CWIPI library" FORCE)
        SET(CWIPI_LIBRARY_FVMC fvmc CACHE FILEPATH
            "CWIPI library" FORCE)
        SET(CWIPI_LIBRARY_BFTC bftc CACHE FILEPATH
            "CWIPI library" FORCE)

        SET(CWIPI_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "CWIPI include" FORCE)

        LINK_DIRECTORIES(${TPDIST}/lib)

        MESSAGE(STATUS "Build CWIPI: ${TPDIST}/lib/lib${CWIPI_LIBRARY}.so ${TPDIST}/lib/lib${CWIPI_LIBRARY_FVMC}.so ${TPDIST}/lib/lib${CWIPI_LIBRARY_BFTC}.so")
        SET(CWIPI_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE ()
        ADD_CUSTOM_TARGET(cwipi-0.8.2 ALL)
        MESSAGE(STATUS "Found CWIPI: ${CWIPI_LIBRARY} ${CWIPI_LIBRARY_FVMC} ${CWIPI_LIBRARY_BFTC}")
        SET(CWIPI_CONFIG_INCLUDE_DIR ${CWIPI_INCLUDE_DIR})
    ENDIF()
ENDIF( NEKTAR_USE_CWIPI )

INCLUDE_DIRECTORIES(SYSTEM ${CWIPI_INCLUDE_DIR})

MARK_AS_ADVANCED(CWIPI_LIBRARY)
MARK_AS_ADVANCED(CWIPI_LIBRARY_FVMC)
MARK_AS_ADVANCED(CWIPI_LIBRARY_BFTC)
MARK_AS_ADVANCED(CWIPI_INCLUDE_DIR)
