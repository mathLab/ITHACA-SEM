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
            cwipi-0.8.2
            URL ${TPURL}/cwipi-0.8.2.tgz
            URL_MD5 "cd28dbea20a08d71f5ff4b4770867268"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/cwipi-0.8.2
            PATCH_COMMAND patch -p 0 < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/cwipi_fix-compile.patch
            BINARY_DIR ${TPBUILD}/cwipi-0.8.2
            TMP_DIR ${TPBUILD}/cwipi-0.8.2-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND
                OMPI_FC=${CMAKE_Fortran_COMPILER}
                OMPI_CC=${CMAKE_C_COMPILER}
                OMPI_CXX=${CMAKE_CXX_COMPILER}
                CFLAGS=-std=c99
                CXXFLAGS=-std=c++11
                ${TPSRC}/cwipi-0.8.2/configure
                CC=${MPI_C_COMPILER}
                CXX=${MPI_CXX_COMPILER}
                FC=${MPI_Fortran_COMPILER}
                --prefix=${TPDIST}
                --libdir=${TPDIST}/lib
                --quiet
            BUILD_COMMAND make -j 1
        )

        THIRDPARTY_LIBRARY(CWIPI_LIBRARY SHARED cwipi
            DESCRIPTION "CWIPI main library")
        THIRDPARTY_LIBRARY(CWIPI_LIBRARY_FVMC SHARED fvmc
            DESCRIPTION "CWIPI fvmc library")
        THIRDPARTY_LIBRARY(CWIPI_LIBRARY_BFTC SHARED bftc
            DESCRIPTION "CWIPI bftc library")

        SET(CWIPI_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "CWIPI include" FORCE)

        LINK_DIRECTORIES(${TPDIST}/lib)

        MESSAGE(STATUS "Build CWIPI:
            ${CWIPI_LIBRARY}
            ${CWIPI_LIBRARY_FVMC}
            ${CWIPI_LIBRARY_BFTC}")
        SET(CWIPI_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE ()
        ADD_CUSTOM_TARGET(cwipi-0.8.2 ALL)
        MESSAGE(STATUS "Found CWIPI:
            ${CWIPI_LIBRARY}
            ${CWIPI_LIBRARY_FVMC}
            ${CWIPI_LIBRARY_BFTC}")
        SET(CWIPI_CONFIG_INCLUDE_DIR ${CWIPI_INCLUDE_DIR})
    ENDIF()
ENDIF( NEKTAR_USE_CWIPI )

INCLUDE_DIRECTORIES(SYSTEM ${CWIPI_INCLUDE_DIR})

MARK_AS_ADVANCED(CWIPI_LIBRARY)
MARK_AS_ADVANCED(CWIPI_LIBRARY_FVMC)
MARK_AS_ADVANCED(CWIPI_LIBRARY_BFTC)
MARK_AS_ADVANCED(CWIPI_INCLUDE_DIR)
