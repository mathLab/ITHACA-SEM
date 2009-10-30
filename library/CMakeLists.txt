CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
PROJECT(Nektar++)

INCLUDE(CMakeExportBuildSettings)

# Use the CMake 2.4 policy for CMake 2.6
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)

# Find the modules included with Nektar
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})
SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/dist CACHE FILEPATH "" FORCE)

SET(NEKTAR_USE_METIS ON CACHE BOOL "Use Metis routines for partitioning and reordering graphs.")
MARK_AS_ADVANCED(NEKTAR_USE_METIS)
SET(NEKTAR_USE_TINYXML_STL ON CACHE BOOL "Use STL with TinyXML library.")
SET(NEKTAR_USE_DIRECT_BLAS_CALLS OFF CACHE BOOL "Use direct BLAS calls rather than NekMatrix calls in the StdRegions and LocalRegions sublibrary")
SET(NEKTAR_FULL_DEBUG OFF CACHE BOOL "Enable Full Debugging.")
SET(NEKTAR_USE_PRECOMPILED_HEADERS OFF CACHE BOOL "Enable the use of precompiled headers.")
#SET(CMAKE_BUILD_TYPE Release CACHE String "")

IF (${CMAKE_COMPILER_IS_GNUCXX}) 
    SET(NEKTAR_ENABLE_PROFILE OFF CACHE BOOL "Compiles Nektar++ with the -pg compiler flag")
    MARK_AS_ADVANCED(NEKTAR_ENABLE_PROFILE)
ENDIF (${CMAKE_COMPILER_IS_GNUCXX}) 

SET(NEKTAR_BUILD_DEMOS ON CACHE BOOL "Build Demo Projects")
SET(NEKTAR_BUILD_UNIT_TESTS OFF CACHE BOOL "Build Unit Tests Projects")

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008" )
	SET(ExpressionTemplateDefault OFF)
ELSE( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008" )
	SET(ExpressionTemplateDefault ON)
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008" )

SET(NEKTAR_USE_EXPRESSION_TEMPLATES ${ExpressionTemplateDefault} CACHE BOOL "")

MARK_AS_ADVANCED(NEKTAR_USE_DIRECT_BLAS_CALLS)
MARK_AS_ADVANCED(NEKTAR_FULL_DEBUG)
MARK_AS_ADVANCED(NEKTAR_USE_MEMORY_POOLS)

MARK_AS_ADVANCED(NEKTAR_USE_PRECOMPILED_HEADERS)
MARK_AS_ADVANCED(NEKTAR_USE_EXPRESSION_TEMPLATES)

SET(NEKTAR_USE_ACML OFF CACHE BOOL "Use the AMD Core Math Library (ACML) for BLAS and Lapack support." )
#SET(NEKTAR_USE_ATLAS OFF CACHE BOOL "Use ATLAS for BLAS and Lapack support." )
SET(NEKTAR_USE_MEMORY_POOLS ON CACHE BOOL "Use memory pools to accelerate memory allocation.")
SET(NEKTAR_USE_BLAS_LAPACK ON CACHE BOOL "Use Blas and lapack routines.")
MARK_AS_ADVANCED(NEKTAR_USE_BLAS_LAPACK)
SET(NEKTAR_USE_NIST_SPARSE_BLAS_TOOLKIT ON CACHE BOOL "Use NIST Fortran Sparse Blas Toolkit library" )
MARK_AS_ADVANCED(NEKTAR_USE_NIST_SPARSE_BLAS_TOOLKIT)
SET(NEKTAR_USE_MKL OFF CACHE BOOL "Use the Intel Math Kernel Library (MKL) for BLAS and Lapack support." )
IF( APPLE )
    SET(NEKTAR_USE_FRAMEWORK ON CACHE BOOL "Use the Mac Framework system for BLAS and Lapack support." )
ENDIF(APPLE)
IF( UNIX )
  SET(NEKTAR_USE_SYSTEM_BLAS_LAPACK OFF CACHE BOOL "Use the system provided blas and lapack libraries")
ENDIF(UNIX)

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

INCLUDE(NektarCommon)


INCLUDE (FindLoki)
INCLUDE (FindTinyXml)

IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
	SET(Boost_USE_STATIC_LIBS OFF)
ELSE(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
	SET(Boost_USE_STATIC_LIBS ON)
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")

SET(Boost_USE_MULTITHREAD ON)
SET(BOOST_ROOT ${CMAKE_SOURCE_DIR}/../ThirdParty/boost 
  ${CMAKE_SOURCE_DIR}/../../ThirdParty/boost)
SET(Boost_ADDITIONAL_VERSIONS "1.35.0" "1.35")
FIND_PACKAGE( Boost 1.35.0 REQUIRED COMPONENTS thread iostreams zlib date_time)


INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${LOKI_INCLUDE_DIR} ${TINYXML_BASE})



IF( NEKTAR_USE_PRECOMPILED_HEADERS )
    ADD_DEFINITIONS(-DNEKTAR_USE_PRECOMPILED_HEADERS)
ELSE( NEKTAR_USE_PRECOMPILED_HEADERS )
    REMOVE_DEFINITIONS(-DNEKTAR_USE_PRECOMPILED_HEADERS)
ENDIF( NEKTAR_USE_PRECOMPILED_HEADERS )





SET(LibrarySubDirs LibUtilities LocalRegions MultiRegions SpatialDomains
StdRegions Stochastic)
SET(UnitTestSubDirs UnitTests)
SET(DemoSubDirs Demos)

SUBDIRS(${LibrarySubDirs} )

IF( NEKTAR_BUILD_UNIT_TESTS )
    SUBDIRS(${UnitTestSubDirs})
ENDIF( NEKTAR_BUILD_UNIT_TESTS )

IF( NEKTAR_BUILD_DEMOS )
    SUBDIRS(${DemoSubDirs})
ENDIF( NEKTAR_BUILD_DEMOS )

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/Nektar++Config.cmake.in 
               ${CMAKE_BINARY_DIR}/Nektar++Config.cmake @ONLY)
INSTALL(FILES ${CMAKE_BINARY_DIR}/Nektar++Config.cmake 
              ${CMAKE_SOURCE_DIR}/FindBoost.cmake
              ${CMAKE_SOURCE_DIR}/FindLoki.cmake
	      ${CMAKE_SOURCE_DIR}/FindFramework.cmake
              ${CMAKE_SOURCE_DIR}/FindACML.cmake
              #${CMAKE_SOURCE_DIR}/FindAtlas.cmake
              ${CMAKE_SOURCE_DIR}/FindNistSparseBlas.cmake
              ${CMAKE_SOURCE_DIR}/FindNativeBlasLapack.cmake
              ${CMAKE_SOURCE_DIR}/FindMKL.cmake
              ${CMAKE_SOURCE_DIR}/FindMetis.cmake
              ${CMAKE_SOURCE_DIR}/FindWin32Lapack.cmake
	      ${CMAKE_SOURCE_DIR}/NektarCommon.cmake
              ${CMAKE_SOURCE_DIR}/FindTinyXml.cmake
              DESTINATION ${CMAKE_INSTALL_PREFIX}) 

INSTALL(EXPORT Nektar++Libraries DESTINATION lib)

CMAKE_EXPORT_BUILD_SETTINGS("Nektar++BuildSettings.cmake")

