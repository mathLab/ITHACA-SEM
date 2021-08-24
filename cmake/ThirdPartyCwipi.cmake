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
    find_package(Cwipi)
    include_directories(${CWIPI_INCLUDE_DIRS})
    ADD_DEFINITIONS(-DNEKTAR_USE_CWIPI)

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

        IF(NOT CMAKE_Fortran_COMPILER)
            MESSAGE(ERROR "MPI_Fortran_COMPILER not set")
        ENDIF()

        IF(NOT NEKTAR_USE_MPI)
            MESSAGE(ERROR "NEKTAR_USE_MPI not set")
        ENDIF()

        EXTERNALPROJECT_ADD(
            cwipi-0.11.1
            URL ${TPURL}/cwipi-0.11.1.tgz
            URL_MD5 "bf51abe03a1007fd084cea670d79d2d5"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/cwipi-0.11.1
            BINARY_DIR ${TPBUILD}/cwipi-0.11.1
            TMP_DIR ${TPBUILD}/cwipi-0.11.1-tmp
            INSTALL_DIR ${TPDIST}
            PATCH_COMMAND patch -p1 < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/cwipi-disable-warnings.patch
            COMMAND patch -p1 < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/cwipi-disable-fortran.patch
            CONFIGURE_COMMAND
                CFLAGS=-w
                CXXFLAGS=-w
                CC=${MPI_C_COMPILER}
                CXX=${MPI_CXX_COMPILER}
                FC=${MPI_Fortran_COMPILER}
                ${CMAKE_COMMAND}
                    -DCMAKE_INSTALL_PREFIX=${TPDIST}
                    ${TPSRC}/cwipi-0.11.1
        )

        THIRDPARTY_LIBRARY(CWIPI_LIBRARY SHARED cwp
            DESCRIPTION "CWIPI main library")

        SET(CWIPI_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "CWIPI include" FORCE)

        MESSAGE(STATUS "Build CWIPI: ${CWIPI_LIBRARY}")
        SET(CWIPI_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE ()
        ADD_CUSTOM_TARGET(cwipi-0.11.1 ALL)
        MESSAGE(STATUS "Found CWIPI: ${CWIPI_LIBRARY}")
        SET(CWIPI_CONFIG_INCLUDE_DIR ${CWIPI_INCLUDE_DIR})
    ENDIF()
ENDIF( NEKTAR_USE_CWIPI )

INCLUDE_DIRECTORIES(SYSTEM ${CWIPI_INCLUDE_DIR})

MARK_AS_ADVANCED(CWIPI_LIBRARY)
MARK_AS_ADVANCED(CWIPI_INCLUDE_DIR)
