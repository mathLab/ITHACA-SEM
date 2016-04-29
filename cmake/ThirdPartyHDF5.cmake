########################################################################
#
# ThirdParty configuration for Nektar++
#
# HDF5
#
########################################################################

OPTION(NEKTAR_USE_HDF5
    "Enable HDF5 I/O support." OFF)

IF (NEKTAR_USE_HDF5)
    INCLUDE(FindHDF5)

    IF (NOT NEKTAR_USE_MPI)
        MESSAGE(FATAL_ERROR "HDF5 requires Nektar++ to be compiled using NEKTAR_USE_MPI.")
    ENDIF()

    IF (HDF5_FOUND)
        IF (NOT HDF5_IS_PARALLEL)
            MESSAGE(FATAL_ERROR "HDF5 detected but is not compiled in parallel.")
        ENDIF()

        INCLUDE_DIRECTORIES(SYSTEM ${HDF5_INCLUDE_DIRS})
        SET(HDF5_CONFIG_INCLUDE_DIR ${TPINC})
        ADD_CUSTOM_TARGET(hdf5-1.8.16 ALL)
    ELSE()
        IF (NOT CMAKE_VERSION VERSION_GREATER 3.1.0)
            MESSAGE(FATAL_ERROR "HDF5 compilation requires CMake 3.1.0 or later.")
        ENDIF()

        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            hdf5-1.8.16
            PREFIX ${TPSRC}
            URL http://ae-nektar.ae.ic.ac.uk/~dmoxey/hdf5-1.8.16.tar.bz2
            URL_MD5 79c1593573ebddf734eee8d43ecfe483
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/hdf5-1.8.16
            BINARY_DIR ${TPBUILD}/hdf5-1.8.16
            TMP_DIR ${TPBUILD}/hdf5-1.8.16-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                -DHDF5_ENABLE_PARALLEL=ON
                -DHDF5_BUILD_CPP_LIB=OFF
                -DBUILD_TESTING=OFF
                -DHDF5_BUILD_TOOLS=OFF
                ${TPSRC}/hdf5-1.8.16
            )

        SET(HDF5_LIBRARIES hdf5-shared CACHE FILEPATH
            "HDF5 libraries" FORCE)
        SET(HDF5_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "HDF5 include directory" FORCE)
        CONSTRUCT_LIBNAME(HDF5_LIBRARIES)

        MARK_AS_ADVANCED(HDF5_LIBRARIES)
        MARK_AS_ADVANCED(HDF5_INCLUDE_DIR)

        LINK_DIRECTORIES(${TPDIST}/lib)

        MESSAGE(STATUS "Build HDF5: ${HDF5_LIBRARIES}")

        SET(HDF5_CONFIG_INCLUDE_DIR ${TPINC})
    ENDIF()
ENDIF (NEKTAR_USE_HDF5)