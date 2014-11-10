# - Nektar++ Config File
#
# Use this module by invoking find_package with the form:
#  FIND_PACKAGE(Nektar++
#    [version] [EXACT]      # Minimum or EXACT version e.g. 1.36.0
#    [REQUIRED]             # Fail with error if Nektar++ is not found
#    )                      #
#
# This sets the following variables:
#  Nektar++_FOUND             - True if headers and requested libraries were found
#  Nektar++_VERSION           - Nektar++_VERSION
#  Nektar++_INCLUDE_DIRS      - Nektar++ include directories
#  Nektar++_LIBRARY_DIRS      - Link directories for Nektar++ libraries
#  Nektar++_DEFINITIONS       - Nektar++ build flags
#  Nektar++_LIBRARIES         - Nektar++ component libraries to be linked
#
#  Nektar++_TP_INCLUDE_DIRS   - Nektar++ ThirdParty include directories
#  Nektar++_TP_LIBRARY_DIRS   - Link directories for Nektar++ ThirdParty libraries
#  Nektar++_TP_LIBRARIES      - Nektar++ ThirdParty libraries to be linked
#
# Example Use:
#  FIND_PACKAGE(Nektar++ REQUIRED)
#  ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})
#  INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
#  LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS} ${NEKTAR++_TP_LIBRARY_DIRS})
#  TARGET_LINK_LIBRARIES(${ProjectName} ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
#

# set basic variables
SET(NEKTAR++_FOUND "ON")
SET(NEKTAR++_VERSION "3.4.0")
SET(NEKTAR++_ROOT_DIR "/tmp/ybao/nektar++/builds/dist")
SET(NEKTAR++_INCLUDE_DIRS "${NEKTAR++_ROOT_DIR}/include/nektar++")
SET(NEKTAR++_LIBRARY_DIRS "${NEKTAR++_ROOT_DIR}/lib64")
SET(NEKTAR++_DEFINITIONS "")
SET(NEKTAR++_LIBRARIES SolverUtils MultiRegions LocalRegions SpatialDomains StdRegions LibUtilities)

SET(NEKTAR++_TP_INCLUDE_DIRS "")
SET(NEKTAR++_TP_LIBRARIES "")
SET(NEKTAR++_TP_LIBRARY_DIRS "/tmp/ybao/nektar++/ThirdParty/dist/lib/")

# add Nektar++_ROOT_DIR to the cmake search path, so Nektars custom FindXXX modules can be used
SET(CMAKE_MODULE_PATH ${NEKTAR++_ROOT_DIR} ${CMAKE_MODULE_PATH})

# set nektars config variables
SET(NEKTAR_USE_MEMORY_POOLS "ON")
IF( NEKTAR_USE_MEMORY_POOLS )
    SET(NEKTAR++_DEFINITIONS "${NEKTAR++_DEFINITIONS} -DNEKTAR_MEMORY_POOL_ENABLED")
ENDIF( NEKTAR_USE_MEMORY_POOLS )

SET(NEKTAR_USE_SMV "OFF")
IF( NEKTAR_USE_SMV )
    SET(NEKTAR++_DEFINITIONS "${NEKTAR++_DEFINITIONS} -DNEKTAR_USING_SMV")
ENDIF( NEKTAR_USE_SMV )

SET(NEKTAR_USE_ACML "OFF")

SET(NEKTAR_USE_ACML "OFF")
IF( NEKTAR_USE_EXPRESSION_TEMPLATES )
    SET(NEKTAR++_DEFINITIONS "${NEKTAR++_DEFINITIONS} -DNEKTAR_USE_EXPRESSION_TEMPLATES")
ENDIF( NEKTAR_USE_EXPRESSION_TEMPLATES )

SET(NEKTAR_USE_WIN32_LAPACK "OFF")
IF( NEKTAR_USE_WIN32_LAPACK )
    GET_FILENAME_COMPONENT(WIN32_BLAS_PATH "" PATH)
    GET_FILENAME_COMPONENT(WIN32_LAPACK_PATH "" PATH)
ENDIF( NEKTAR_USE_WIN32_LAPACK )

SET(NEKTAR_USE_BLAS_LAPACK "ON")
IF( NEKTAR_USE_BLAS_LAPACK )
    SET(NEKTAR++_DEFINITIONS "${NEKTAR++_DEFINITIONS} -DNEKTAR_USING_LAPACK -DNEKTAR_USING_BLAS")
ENDIF( NEKTAR_USE_BLAS_LAPACK )

SET(NEKTAR_FULL_DEBUG "OFF")
IF ( NEKTAR_FULL_DEBUG )
    SET(NEKTAR++_DEFINITIONS "${NEKTAR++_DEFINITIONS} -DNEKTAR_FULLDEBUG")
ENDIF( NEKTAR_FULL_DEBUG)

SET(Boost_INCLUDE_DIRS "/apps/boost/1.49.0/include")
SET(Boost_LIBRARIES "/apps/boost/1.49.0/lib/libboost_thread.so;/apps/boost/1.49.0/lib/libboost_iostreams.so;/apps/boost/1.49.0/lib/libboost_date_time.so;/apps/boost/1.49.0/lib/libboost_filesystem.so;/apps/boost/1.49.0/lib/libboost_system.so;/apps/boost/1.49.0/lib/libboost_program_options.so;/apps/boost/1.49.0/lib/libboost_regex.so")
SET(Boost_LIBRARY_DIRS "/apps/boost/1.49.0/lib")
SET(NEKTAR++_TP_INCLUDE_DIRS ${NEKTAR++_TP_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
SET(NEKTAR++_TP_LIBRARIES ${NEKTAR++_TP_LIBRARIES} ${Boost_LIBRARIES})

SET(NEKTAR_USE_MPI "ON")
SET(MPI_LIBRARY "/usr/lib/openmpi/lib/libmpi_cxx.so")
SET(MPI_EXTRA_LIBRARY "/usr/lib/openmpi/lib/libmpi.so;/usr/lib/openmpi/lib/libopen-rte.so;/usr/lib/openmpi/lib/libopen-pal.so;/usr/lib/x86_64-linux-gnu/libdl.so;/usr/lib/x86_64-linux-gnu/libnsl.so;/usr/lib/x86_64-linux-gnu/libutil.so;/usr/lib/x86_64-linux-gnu/libm.so;/usr/lib/x86_64-linux-gnu/libdl.so")
SET(MPI_INCLUDE_PATH "/usr/lib/openmpi/include;/usr/lib/openmpi/include/openmpi")
SET(NEKTAR++_TP_INCLUDE_DIRS ${NEKTAR++_TP_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
SET(NEKTAR++_TP_LIBRARIES ${NEKTAR++_TP_LIBRARIES} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
IF( NEKTAR_USE_MPI )
    SET(NEKTAR++_DEFINITIONS "${NEKTAR++_DEFINITIONS} -DNEKTAR_USE_MPI")
ENDIF( NEKTAR_USE_MPI )

SET(LOKI_ADDITIONAL_INCLUDE_DIRS "/tmp/ybao/nektar++/ThirdParty/loki-0.1.3/include")
SET(NEKTAR++_TP_INCLUDE_DIRS ${NEKTAR++_TP_INCLUDE_DIRS} ${LOKI_ADDITIONAL_INCLUDE_DIRS})

SET(TINYXML_ADDITIONAL_INCLUDE_DIRS "")
SET(NEKTAR++_TP_INCLUDE_DIRS ${NEKTAR++_TP_INCLUDE_DIRS} ${TINYXML_ADDITIONAL_INCLUDE_DIRS})

SET(FFTW_INCLUDE_DIR "/apps/fftw/3.3.2/lib/../include")
SET(NEKTAR++_TP_INCLUDE_DIRS ${NEKTAR++_TP_INCLUDE_DIRS} ${FFTW_INCLUDE_DIR})

SET(ARPACK_INCLUDE_DIR "")
SET(NEKTAR++_TP_INCLUDE_DIRS ${NEKTAR++_TP_INCLUDE_DIRS} ${ARPACK_INCLUDE_DIR})

SET(VTK_INCLUDE_DIRS "")
SET(NEKTAR++_TP_INCLUDE_DIRS ${NEKTAR++_TP_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})

# find and add Nektar++ libraries
INCLUDE(${NEKTAR++_LIBRARY_DIRS}/Nektar++Libraries.cmake)

# platform dependent options
if(${CMAKE_SYSTEM} MATCHES "Linux.*")
    set(NEKTAR++_TP_LIBRARIES ${NEKTAR++_TP_LIBRARIES} rt)
    SET(NEKTAR++_DEFINITIONS "${NEKTAR++_DEFINITIONS} -pthread")
endif(${CMAKE_SYSTEM} MATCHES "Linux.*")

if(${CMAKE_SYSTEM} MATCHES "Darwin-*")
    SET(NEKTAR++_DEFINITIONS "${NEKTAR++_DEFINITIONS} -Wl,-undefined,dynamic_lookup -Wl,-rpath,${CMAKE_INSTALL_PREFIX}/${LIB_DIR} -Wl,-rpath,${Boost_LIBRARY_DIRS}")
endif(${CMAKE_SYSTEM} MATCHES "Darwin-*")
