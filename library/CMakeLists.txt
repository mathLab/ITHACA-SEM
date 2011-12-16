CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
PROJECT(Nektar++)

INCLUDE(CMakeExportBuildSettings)
ENABLE_TESTING()

# Use the CMake 2.4 policy for CMake 2.6
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)

# Find the modules included with Nektar
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})
IF( NOT NEKTAR_PASSED_FIRST_CONFIGURE )
  SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/dist CACHE PATH "" FORCE)
ENDIF()

SET(LIB_DIR "lib")
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
	SET(LIB_DIR "lib64")
ENDIF()

# Set up RPATH
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_DIR}")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${LIB_DIR}")
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
	"${CMAKE_INSTALL_PREFIX}/${LIB_DIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
	SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_DIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")


SET(NEKTAR_USE_METIS ON CACHE BOOL "Use Metis routines for partitioning and reordering graphs.")
MARK_AS_ADVANCED(NEKTAR_USE_METIS)
SET(NEKTAR_USE_ARPACK OFF CACHE BOOL "Use Arpack routines for evaluating the
eigenvalues and eigenvectors")
MARK_AS_ADVANCED(NEKTAR_USE_ARPACK)

SET(NEKTAR_USE_FFTW OFF CACHE BOOL "Use FFTW routines for performing the Fast Fourier Transform.")
MARK_AS_ADVANCED(NEKTAR_USE_FFTW)


# Parallelisation
SET(NEKTAR_USE_MPI OFF CACHE BOOL "Use MPICH2 for parallelisation.")
SET(MPI_SYNC OFF CACHE BOOL "Use synchronous MPI calls when sending data.")
MARK_AS_ADVANCED(MPI_SYNC)
IF( NEKTAR_USE_MPI )
    INCLUDE (FindGSLib)
ENDIF( NEKTAR_USE_MPI )

SET(NEKTAR_USE_TINYXML_STL ON CACHE BOOL "Use STL with TinyXML library.")
SET(NEKTAR_USE_DIRECT_BLAS_CALLS OFF CACHE BOOL "Use direct BLAS calls rather than NekMatrix calls in the StdRegions and LocalRegions sublibrary")
SET(NEKTAR_FULL_DEBUG OFF CACHE BOOL "Enable Full Debugging.")
SET(NEKTAR_USE_PRECOMPILED_HEADERS OFF CACHE BOOL "Enable the use of precompiled headers.")

IF (${CMAKE_COMPILER_IS_GNUCXX}) 
    SET(NEKTAR_ENABLE_PROFILE OFF CACHE BOOL "Compiles Nektar++ with the -pg compiler flag")
    MARK_AS_ADVANCED(NEKTAR_ENABLE_PROFILE)
ENDIF (${CMAKE_COMPILER_IS_GNUCXX}) 

SET(NEKTAR_BUILD_DEMOS ON CACHE BOOL "Build Demo Projects")
SET(NEKTAR_BUILD_UNIT_TESTS OFF CACHE BOOL "Build Unit Tests Projects")
SET(NEKTAR_BUILD_SHARED ON CACHE BOOL "Build Nektar++ as shared libraries")
SET(NEKTAR_BUILD_TIMINGS OFF CACHE BOOL "Build Timings directory")
MARK_AS_ADVANCED(NEKTAR_BUILD_TIMINGS)

SET(ExpressionTemplateDefault OFF)

SET(NEKTAR_USE_EXPRESSION_TEMPLATES ${ExpressionTemplateDefault} CACHE BOOL "")

MARK_AS_ADVANCED(NEKTAR_USE_DIRECT_BLAS_CALLS)
MARK_AS_ADVANCED(NEKTAR_FULL_DEBUG)
MARK_AS_ADVANCED(NEKTAR_USE_MEMORY_POOLS)

MARK_AS_ADVANCED(NEKTAR_USE_PRECOMPILED_HEADERS)
MARK_AS_ADVANCED(NEKTAR_USE_EXPRESSION_TEMPLATES)

SET(NEKTAR_USE_ACML OFF CACHE BOOL "Use the AMD Core Math Library (ACML) for BLAS and Lapack support." )
SET(NEKTAR_USE_MEMORY_POOLS ON CACHE BOOL "Use memory pools to accelerate memory allocation.")
SET(NEKTAR_USE_BLAS_LAPACK ON CACHE BOOL "Use Blas and lapack routines.")
MARK_AS_ADVANCED(NEKTAR_USE_BLAS_LAPACK)
SET(NEKTAR_USE_NIST_SPARSE_BLAS_TOOLKIT ON CACHE BOOL "Use NIST Fortran Sparse Blas Toolkit library" )
MARK_AS_ADVANCED(NEKTAR_USE_NIST_SPARSE_BLAS_TOOLKIT)
SET(NEKTAR_USE_MKL OFF CACHE BOOL "Use the Intel Math Kernel Library (MKL) for BLAS and Lapack support." )
IF( APPLE )
    SET(NEKTAR_USE_ACCELERATE_FRAMEWORK ON CACHE BOOL "Use the Mac Accelerate Framework for BLAS and Lapack support." )
    SET(NEKTAR_USE_CHUD_FRAMEWORK OFF CACHE BOOL "Use the Mac CHUD Framework (for Shark support)." )
    MARK_AS_ADVANCED(NEKTAR_USE_CHUD_FRAMEWORK)
ELSE(APPLE)
        IF( UNIX )
          SET(NEKTAR_USE_SYSTEM_BLAS_LAPACK ON CACHE BOOL "Use the system provided blas and lapack libraries")
          ENDIF(UNIX)
ENDIF(APPLE)

# Turn on NEKTAR_USE_WIN32_LAPACK if we are in Windows and the libraries exist.
IF( WIN32 )
	FIND_LIBRARY( WIN32_BLAS NAMES blas_win32 PATHS ${CMAKE_SOURCE_DIR}/../ThirdParty )
	FIND_LIBRARY( WIN32_LAPACK NAMES lapack_win32 PATHS ${CMAKE_SOURCE_DIR}/../ThirdParty )
	IF( WIN32_BLAS AND WIN32_LAPACK )
		SET(Win32LapackDefault ON)
    ELSE(WIN32_BLAS AND WIN32_LAPACK)
		SET(Win32LapackDefault OFF)
    ENDIF(WIN32_BLAS AND WIN32_LAPACK)
ELSE(WIN32)
	SET(Win32LapackDefault OFF)
ENDIF(WIN32)

SET(NEKTAR_USE_WIN32_LAPACK ${Win32LapackDefault} CACHE BOOL "Use Win32 Lapack provided with the Third Party Distribution.")

INCLUDE (NektarCommon)

INCLUDE (FindLoki)
INCLUDE (FindTinyXml)

IF(NEKTAR_USE_FFTW) 
    INCLUDE (FindFFTW)
ENDIF(NEKTAR_USE_FFTW)

IF(NEKTAR_USE_ARPACK)
    INCLUDE (FindArpack)
ENDIF(NEKTAR_USE_ARPACK)

IF( NOT WIN32 )
    SET(Boost_USE_STATIC_LIBS OFF)
ELSE( NOT WIN32 )
    SET(Boost_USE_STATIC_LIBS ON)
ENDIF( NOT WIN32 )

SET(Boost_USE_MULTITHREAD ON)
SET(Boost_ADDITIONAL_VERSIONS "1.48.0" "1.47.0" "1.47" "1.46" "1.46.1" "1.40" "1.40.0" "1.35.0" "1.35")

IF( NOT BOOST_ROOT )
    #If the user has not set BOOST_ROOT, look in a couple common places first.
    SET(BOOST_ROOT ${CMAKE_SOURCE_DIR}/../ThirdParty/boost) 
    FIND_PACKAGE( Boost 1.40.0 COMPONENTS thread iostreams zlib date_time filesystem system)
    SET(BOOST_ROOT ${CMAKE_SOURCE_DIR}/../../ThirdParty/boost)
    FIND_PACKAGE( Boost 1.40.0 COMPONENTS thread iostreams zlib date_time filesystem system)
ELSE()
    FIND_PACKAGE( Boost 1.48.0 COMPONENTS thread iostreams zlib date_time filesystem system)
ENDIF()

IF(NOT Boost_ZLIB_FOUND)
	FIND_PACKAGE( ZLIB )
	IF (ZLIB_FOUND)
		SET(Boost_ZLIB_LIBRARY ${ZLIB_LIBRARIES})
		SET(Boost_ZLIB_LIBRARY_RELEASE ${ZLIB_LIBRARIES})
		SET(Boost_ZLIB_LIBRARY_DEBUG ${ZLIB_LIBRARIES})
	ENDIF()
ENDIF()

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${LOKI_INCLUDE_DIR} ${TINYXML_BASE})

IF(NEKTAR_USE_FFTW) 
    INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
ENDIF(NEKTAR_USE_FFTW)

IF( NEKTAR_USE_MPI )
    INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )
ENDIF( NEKTAR_USE_MPI )

IF(NEKTAR_USE_ARPACK)
    INCLUDE_DIRECTORIES(${ARPACK_INCLUDE_DIR})
ENDIF(NEKTAR_USE_ARPACK)

IF( NEKTAR_USE_PRECOMPILED_HEADERS )
    ADD_DEFINITIONS(-DNEKTAR_USE_PRECOMPILED_HEADERS)
ELSE( NEKTAR_USE_PRECOMPILED_HEADERS )
    REMOVE_DEFINITIONS(-DNEKTAR_USE_PRECOMPILED_HEADERS)
ENDIF( NEKTAR_USE_PRECOMPILED_HEADERS )


IF( NEKTAR_USE_NIST_SPARSE_BLAS_TOOLKIT )
    INCLUDE (FindNistSparseBlas)
ENDIF( NEKTAR_USE_NIST_SPARSE_BLAS_TOOLKIT )


SET(LibrarySubDirs LibUtilities LocalRegions MultiRegions SpatialDomains
StdRegions)
SET(UnitTestSubDirs UnitTests)
SET(DemoSubDirs Demos)
SET(TimingsSubDirs Timings)

SUBDIRS(${LibrarySubDirs} )

IF( NEKTAR_BUILD_UNIT_TESTS )
    SUBDIRS(${UnitTestSubDirs})
ENDIF( NEKTAR_BUILD_UNIT_TESTS )

IF( NEKTAR_BUILD_DEMOS )
    SUBDIRS(${DemoSubDirs})
ENDIF( NEKTAR_BUILD_DEMOS )

IF( NEKTAR_BUILD_TIMINGS ) 
    SUBDIRS(${TimingsSubDirs})
ENDIF( NEKTAR_BUILD_TIMINGS ) 

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/Nektar++Config.cmake.in 
               ${CMAKE_BINARY_DIR}/Nektar++Config.cmake @ONLY)
             
INSTALL(FILES ${CMAKE_BINARY_DIR}/Nektar++Config.cmake 
              ${CMAKE_SOURCE_DIR}/FindLoki.cmake
              ${CMAKE_SOURCE_DIR}/FindAccelerateFramework.cmake
              ${CMAKE_SOURCE_DIR}/FindCHUDFramework.cmake
              ${CMAKE_SOURCE_DIR}/FindACML.cmake
              ${CMAKE_SOURCE_DIR}/FindArpack.cmake
              ${CMAKE_SOURCE_DIR}/FindNistSparseBlas.cmake
              ${CMAKE_SOURCE_DIR}/FindNativeBlasLapack.cmake
              ${CMAKE_SOURCE_DIR}/FindMKL.cmake
              ${CMAKE_SOURCE_DIR}/FindMetis.cmake
	      ${CMAKE_SOURCE_DIR}/FindFFTW.cmake
              ${CMAKE_SOURCE_DIR}/FindWin32Lapack.cmake
              ${CMAKE_SOURCE_DIR}/NektarCommon.cmake
              ${CMAKE_SOURCE_DIR}/FindTinyXml.cmake
              ${CMAKE_SOURCE_DIR}/FindGSLib.cmake
              DESTINATION ${CMAKE_INSTALL_PREFIX}) 


INSTALL(EXPORT Nektar++Libraries DESTINATION ${LIB_DIR})
SET(NEKTAR_PASSED_FIRST_CONFIGURE ON CACHE BOOL INTERNAL)


# Packaging
SET(NEKTAR_PACKAGE_GENERATOR "None" CACHE STRING "Support Packaging: RPM or DEB" )
MARK_AS_ADVANCED(NEKTAR_PACKAGE_GENERATOR)

# Define components
SET(CPACK_COMPONENTS_ALL dev lib demos extra-demos)

#SET(CPACK_RPM_PACKAGE_DEBUG ON)
SET(CPACK_COMPONENT_LIB_NAME "nektar++-lib")
SET(CPACK_COMPONENT_LIB_DISPLAY_NAME "nektar++-lib")
SET(CPACK_COMPONENT_LIB_DESCRIPTION "Nektar++ Libraries")

SET(CPACK_COMPONENT_DEV_NAME "nektar++-dev")
SET(CPACK_COMPONENT_DEV_DISPLAY_NAME "nektar++-dev")
SET(CPACK_COMPONENT_DEV_DESCRIPTION "Development files for Nektar++")
SET(CPACK_COMPONENT_DEV_DEPENDS lib)

SET(CPACK_COMPONENT_DEMOS_NAME "nektar++-demos")
SET(CPACK_COMPONENT_DEMOS_DISPLAY_NAME "nektar++-demos")
SET(CPACK_COMPONENT_DEMOS_DESCRIPTION "Framework demonstration binaries")
SET(CPACK_COMPONENT_DEMOS_DEPENDS lib)

SET(CPACK_COMPONENT_EXTRA-DEMOS_NAME "nektar++-extra-demos")
SET(CPACK_COMPONENT_EXTRA-DEMOS_DISPLAY_NAME "nektar++-extra-demos")
SET(CPACK_COMPONENT_EXTRA-DEMOS_DESCRIPTION 
	"Framework extra demonstration binaries")
SET(CPACK_COMPONENT_EXTRA-DEMOS_DEPENDS lib)

INCLUDE(NektarPackageGenerator)
