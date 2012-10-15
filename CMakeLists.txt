CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")

PROJECT(Nektar++)

SET(NEKTAR_VERSION_MAJOR 3)
SET(NEKTAR_VERSION_MINOR 2)
SET(NEKTAR_VERSION_PATCH 0)
SET(NEKTAR_VERSION ${NEKTAR_VERSION_MAJOR}.${NEKTAR_VERSION_MINOR}.${NEKTAR_VERSION_PATCH})

INCLUDE(CMakeExportBuildSettings)
ENABLE_TESTING()

# Use the CMake 2.4 policy for CMake 2.6
IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
ENDIF(COMMAND cmake_policy)

# Find the modules included with Nektar
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
IF( NOT NEKTAR_PASSED_FIRST_CONFIGURE )
  SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/dist CACHE PATH "" FORCE)
ENDIF()

SET(LIB_DIR "lib")
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
	SET(LIB_DIR "lib64")
ENDIF()

SET(NEKTAR_BIN_DIR bin)
SET(NEKTAR_LIB_DIR ${LIB_DIR})
SET(NEKTAR_INCLUDE_DIR include/nektar++)
SET(NEKTAR_SHARE_DIR share/nektar++)
SET(NEKTAR_LIBRARY_TYPE "SHARED")


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


# Components of the Nektar++ package to build
SET(NEKTAR_BUILD_LIBRARY     ON CACHE BOOL
    "Build main Nektar++ libraries.")
SET(NEKTAR_BUILD_DEMOS       ON CACHE BOOL
    "Build demonstration codes.")
SET(NEKTAR_BUILD_SOLVERS     ON CACHE BOOL
    "Build example solvers distributed with Nektar++.")
SET(NEKTAR_BUILD_UTILITIES   ON CACHE BOOL
    "Build utilities for pre-processing and post-processing Nektar++ files.")
SET(NEKTAR_BUILD_UNIT_TESTS  ON CACHE BOOL
    "Build unit tests.")
SET(NEKTAR_BUILD_REGRESSION_TESTS ON CACHE BOOL
    "Build regression tests for checking the correctness of the code.")
SET(NEKTAR_BUILD_TIMINGS     OFF CACHE BOOL 
    "Build benchmark timing codes.")

# Build options
SET(NEKTAR_USE_DIRECT_BLAS_CALLS OFF CACHE BOOL 
    "Use direct BLAS calls rather than NekMatrix calls in the StdRegions and LocalRegions sublibrary")
MARK_AS_ADVANCED(NEKTAR_USE_DIRECT_BLAS_CALLS)
SET(NEKTAR_FULL_DEBUG  OFF CACHE BOOL 
    "Enable Full Debugging.")
MARK_AS_ADVANCED(NEKTAR_FULL_DEBUG)
#SET(NEKTAR_USE_PRECOMPILED_HEADERS OFF CACHE BOOL 
#    "Enable the use of precompiled headers.")
#MARK_AS_ADVANCED(NEKTAR_USE_PRECOMPILED_HEADERS)

IF (${CMAKE_COMPILER_IS_GNUCXX}) 
    SET(NEKTAR_ENABLE_PROFILE OFF CACHE BOOL 
        "Compiles Nektar++ with the -pg compiler flag")
    MARK_AS_ADVANCED(NEKTAR_ENABLE_PROFILE)
ENDIF (${CMAKE_COMPILER_IS_GNUCXX}) 

SET(ExpressionTemplateDefault OFF)
SET(NEKTAR_USE_EXPRESSION_TEMPLATES ${ExpressionTemplateDefault} CACHE BOOL "")
MARK_AS_ADVANCED(NEKTAR_USE_EXPRESSION_TEMPLATES)

SET(NEKTAR_USE_ACML OFF CACHE BOOL 
    "Use the AMD Core Math Library (ACML) for BLAS and Lapack support." )
SET(NEKTAR_USE_MEMORY_POOLS ON CACHE BOOL 
    "Use memory pools to accelerate memory allocation.")
MARK_AS_ADVANCED(NEKTAR_USE_MEMORY_POOLS)
SET(NEKTAR_USE_BLAS_LAPACK ON CACHE BOOL 
    "Use Blas and lapack routines.")
SET(NEKTAR_USE_MKL OFF CACHE BOOL 
    "Use the Intel Math Kernel Library (MKL) for BLAS and Lapack support." )

IF( APPLE )
    SET(NEKTAR_USE_ACCELERATE_FRAMEWORK ON CACHE BOOL 
        "Use the Mac Accelerate Framework for BLAS and Lapack support." )
    SET(NEKTAR_USE_CHUD_FRAMEWORK OFF CACHE BOOL 
        "Use the Mac CHUD Framework (for Shark support)." )
    MARK_AS_ADVANCED(NEKTAR_USE_CHUD_FRAMEWORK)
ELSE(APPLE)
    IF( UNIX )
        SET(NEKTAR_USE_SYSTEM_BLAS_LAPACK ON CACHE BOOL 
            "Use the system provided blas and lapack libraries")
    ENDIF(UNIX)
ENDIF(APPLE)

# Turn on NEKTAR_USE_WIN32_LAPACK if we are in Windows and the libraries exist.
IF( WIN32 )
	FIND_LIBRARY( WIN32_BLAS NAMES blas_win32 PATHS ${CMAKE_SOURCE_DIR}/ThirdParty ${CMAKE_SOURCE_DIR}/../ThirdParty )
	FIND_LIBRARY( WIN32_LAPACK NAMES lapack_win32 PATHS ${CMAKE_SOURCE_DIR}/ThirdParty ${CMAKE_SOURCE_DIR}/../ThirdParty )
	IF( WIN32_BLAS AND WIN32_LAPACK )
		SET(Win32LapackDefault ON)
    ELSE(WIN32_BLAS AND WIN32_LAPACK)
		SET(Win32LapackDefault OFF)
    ENDIF(WIN32_BLAS AND WIN32_LAPACK)
ELSE(WIN32)
	SET(Win32LapackDefault OFF)
ENDIF(WIN32)

IF( MSVC )
	# Needed for M_PI to be visible in visual studio.
	ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
	
	# Removes the warnings about unsafe methods such as strcpy, std::copy,
    # memcmp, etc.
	ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE)
ENDIF( )


SET(NEKTAR_USE_WIN32_LAPACK ${Win32LapackDefault} CACHE BOOL 
    "Use Win32 Lapack provided with the Third Party Distribution.")
IF( NOT WIN32 )
    MARK_AS_ADVANCED(NEKTAR_USE_WIN32_LAPACK)
ENDIF( NOT WIN32 )

# Include Nektar++ common macros
INCLUDE (NektarCommon)

# Find Libraries
SET(TPURL http://www.nektar.info/thirdparty)
SET(TPSRC ${CMAKE_SOURCE_DIR}/ThirdParty)
IF (NOT EXISTS ${TPSRC})
    FILE(MAKE_DIRECTORY ${TPSRC})
ENDIF ()
INCLUDE (ThirdPartyTinyxml)
INCLUDE (ThirdPartyLoki)
INCLUDE (ThirdPartyMetis)
INCLUDE (ThirdPartySpblastk)
INCLUDE (ThirdPartyBoost)
INCLUDE (ThirdPartyZlib)
INCLUDE (ThirdPartyFFTW)
INCLUDE (ThirdPartyArpack)
INCLUDE (ThirdPartyMPI)
INCLUDE (ThirdPartyVTK)

IF( NEKTAR_USE_MKL )
    INCLUDE (FindMKL)
    SET(NEKTAR_USING_BLAS TRUE)
    SET(NEKTAR_USING_LAPACK TRUE)
    SET(NEKTAR_USING_MKL TRUE)
ENDIF( NEKTAR_USE_MKL )

IF( NEKTAR_USE_ACCELERATE_FRAMEWORK )
    INCLUDE (FindAccelerateFramework)
ENDIF( NEKTAR_USE_ACCELERATE_FRAMEWORK )

IF( NEKTAR_USE_CHUD_FRAMEWORK )
    INCLUDE (FindCHUDFramework)
    ADD_DEFINITIONS(-DNEKTAR_USING_CHUD)	
ENDIF( NEKTAR_USE_CHUD_FRAMEWORK )

IF( NEKTAR_USE_ACML )
    INCLUDE (FindACML)
    SET(NEKTAR_USING_BLAS TRUE)
    SET(NEKTAR_USING_LAPACK TRUE)
    SET(NEKTAR_USING_ACML TRUE)
    ADD_DEFINITIONS(-DUSE_ACML)
ENDIF( NEKTAR_USE_ACML )

IF( NEKTAR_USE_WIN32_LAPACK )
    INCLUDE (FindWin32Lapack)
ENDIF( NEKTAR_USE_WIN32_LAPACK )

IF( NEKTAR_USE_SYSTEM_BLAS_LAPACK )
    INCLUDE(FindNativeBlasLapack)
    SET(NEKTAR_USING_BLAS TRUE)
    SET(NEKTAR_USING_LAPACK TRUE)
ENDIF( NEKTAR_USE_SYSTEM_BLAS_LAPACK )

IF( NEKTAR_USE_BLAS_LAPACK )
    SET(NEKTAR_USING_BLAS TRUE)
    SET(NEKTAR_USING_LAPACK TRUE)
    ADD_DEFINITIONS( -DNEKTAR_USING_LAPACK -DNEKTAR_USING_BLAS )
ENDIF( NEKTAR_USE_BLAS_LAPACK )

IF( NEKTAR_USE_TINYXML_STL )
    ADD_DEFINITIONS( -DTIXML_USE_STL)
ENDIF( NEKTAR_USE_TINYXML_STL )

IF( NEKTAR_USE_DIRECT_BLAS_CALLS )    
    ADD_DEFINITIONS(-DNEKTAR_USING_DIRECT_BLAS_CALLS)
ENDIF( NEKTAR_USE_DIRECT_BLAS_CALLS )

IF( NEKTAR_USE_EXPRESSION_TEMPLATES )
    ADD_DEFINITIONS(-DNEKTAR_USE_EXPRESSION_TEMPLATES -DNEKTAR_USING_CMAKE)
ENDIF( NEKTAR_USE_EXPRESSION_TEMPLATES )

IF( NEKTAR_USE_MEMORY_POOLS )
    ADD_DEFINITIONS(-DNEKTAR_MEMORY_POOL_ENABLED)
ELSE( NEKTAR_USE_MEMORY_POOLS )
    REMOVE_DEFINITIONS(-DNEKTAR_MEMORY_POOL_ENABLED)
ENDIF( NEKTAR_USE_MEMORY_POOLS )

SET(Boost_USE_STATIC_LIBS OFF)
IF( WIN32 )
    # The auto-linking feature has problems with USE_STATIC_LIBS off, so we use BOOST_ALL_NO_LIB to turn it off.
    # Several boost libraries headers aren't configured correctly if USE_STATIC_LIBS is off, so we explicitly 
    # say they are dynamic with the remaining definitions.
    ADD_DEFINITIONS(-DBOOST_ALL_NO_LIB -DBOOST_PROGRAM_OPTIONS_DYN_LINK -DBOOST_IOSTREAMS_DYN_LINK -DBOOST_THREAD_DYN_LINK)
ENDIF( )


# Set up include directories
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${LOKI_INCLUDE_DIR} ${TINYXML_BASE})

#IF( NEKTAR_USE_PRECOMPILED_HEADERS )
#    ADD_DEFINITIONS(-DNEKTAR_USE_PRECOMPILED_HEADERS)
#ELSE( NEKTAR_USE_PRECOMPILED_HEADERS )
#    REMOVE_DEFINITIONS(-DNEKTAR_USE_PRECOMPILED_HEADERS)
#ENDIF( NEKTAR_USE_PRECOMPILED_HEADERS )

# Build active components
IF (NEKTAR_BUILD_LIBRARY)
    SET(NEKTAR++_LIBRARIES SolverUtils LibUtilities StdRegions SpatialDomains LocalRegions
        MultiRegions)
    INCLUDE_DIRECTORIES(library)
    ADD_SUBDIRECTORY(library)
ENDIF (NEKTAR_BUILD_LIBRARY)

IF (NEKTAR_BUILD_SOLVERS)
    INCLUDE_DIRECTORIES(solvers)
    ADD_SUBDIRECTORY(solvers)
ENDIF (NEKTAR_BUILD_SOLVERS)

IF (NEKTAR_BUILD_UTILITIES)
    INCLUDE_DIRECTORIES(utilities)
    ADD_SUBDIRECTORY(utilities)
ENDIF (NEKTAR_BUILD_UTILITIES)

IF (NEKTAR_BUILD_REGRESSION_TESTS)
    INCLUDE_DIRECTORIES(regressionTests)
    ADD_SUBDIRECTORY(regressionTests)
ENDIF (NEKTAR_BUILD_REGRESSION_TESTS)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/Nektar++Config.cmake.in 
               ${CMAKE_BINARY_DIR}/Nektar++Config.cmake @ONLY)
             
INSTALL(FILES ${CMAKE_BINARY_DIR}/Nektar++Config.cmake 
              ${CMAKE_SOURCE_DIR}/cmake/FindLoki.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindAccelerateFramework.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindCHUDFramework.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindACML.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindArpack.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindNistSparseBlas.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindNativeBlasLapack.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindMKL.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindMetis.cmake
	          ${CMAKE_SOURCE_DIR}/cmake/FindFFTW.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindWin32Lapack.cmake
              ${CMAKE_SOURCE_DIR}/cmake/NektarCommon.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindTinyXml.cmake
              ${CMAKE_SOURCE_DIR}/cmake/FindGSMPI.cmake
              DESTINATION ${CMAKE_INSTALL_PREFIX}) 

INSTALL(EXPORT Nektar++Libraries DESTINATION ${LIB_DIR})

SET(NEKTAR_PASSED_FIRST_CONFIGURE ON CACHE BOOL INTERNAL)
MARK_AS_ADVANCED(NEKTAR_PASSED_FIRST_CONFIGURE)

INCLUDE(NektarPackageGenerator)
