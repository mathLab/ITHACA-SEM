CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
# Toggleable build type
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

PROJECT(Nektar++ C CXX)

CMAKE_POLICY(SET CMP0022 NEW)

# Defer to use _ROOT variables to avoid Zlib warning
IF (POLICY CMP0074)
    CMAKE_POLICY(SET CMP0074 OLD)
ENDIF()

# Nektar++ requires C++11. Try to infer this for older CMake versions (less than
# 3.1.0)
IF ("${CMAKE_VERSION}" VERSION_LESS "3.1")
    IF (NOT MSVC)
        INCLUDE(CheckCXXCompilerFlag)
        CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

        IF (COMPILER_SUPPORTS_CXX11)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        ELSE()
            MESSAGE(FATAL_ERROR "Nektar++ requires a compiler with C++11 support.")
        ENDIF()
    ELSEIF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
        MESSAGE(FATAL_ERROR "Nektar++ requires full C++11 support: please upgrade to Visual Studio 2013 or later")
    ENDIF()
ELSE()
    SET(CMAKE_CXX_STANDARD 11)
    SET(CMAKE_CXX_STANDARD_REQUIRED ON)
    SET(CMAKE_CXX_EXTENSIONS OFF)
ENDIF()

INCLUDE(CheckLanguage)

CHECK_LANGUAGE(Fortran)
IF(CMAKE_Fortran_COMPILER)
    ENABLE_LANGUAGE(Fortran)
ELSE()
    MESSAGE(STATUS "No Fortran support")
ENDIF()

# Helps organize projects in IDEs.
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

# Extract the version number from the VERSION file and set in CMake
# The format of this file must be X.X.X
FILE(STRINGS "VERSION" NEKVER)
STRING(REPLACE "." ";" NEKVERLIST ${NEKVER})
LIST(GET NEKVERLIST 0 NEKTAR_VERSION_MAJOR)
LIST(GET NEKVERLIST 1 NEKTAR_VERSION_MINOR)
LIST(GET NEKVERLIST 2 NEKTAR_VERSION_PATCH)
SET(NEKTAR_VERSION ${NEKTAR_VERSION_MAJOR}.${NEKTAR_VERSION_MINOR}.${NEKTAR_VERSION_PATCH})

# Add support for CMAKE_DEPENDENT_OPTION
INCLUDE(CMakeDependentOption)
INCLUDE(CMakeParseArguments)

# Enable CTest.
ENABLE_TESTING()

# Find the modules included with Nektar
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Set default install path if not provided on the command-line
IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/dist CACHE PATH "" FORCE)
ENDIF()

# Find default search paths for OS X; adapted from Stack Overflow question
# 1487752.
IF (APPLE)
    EXECUTE_PROCESS(COMMAND which port
        RESULT_VARIABLE DETECT_MACPORTS
        OUTPUT_VARIABLE MACPORTS_PREFIX
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    EXECUTE_PROCESS(COMMAND brew --prefix
        RESULT_VARIABLE DETECT_HOMEBREW
        OUTPUT_VARIABLE HOMEBREW_PREFIX
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

    IF (${DETECT_MACPORTS} EQUAL 0)
        GET_FILENAME_COMPONENT(MACPORTS_PREFIX ${MACPORTS_PREFIX} DIRECTORY)
        GET_FILENAME_COMPONENT(MACPORTS_PREFIX ${MACPORTS_PREFIX} DIRECTORY)
        SET(CMAKE_LIBRARY_PATH ${MACPORTS_PREFIX}/lib ${CMAKE_LIBRARY_PATH})
        SET(CMAKE_INCLUDE_PATH ${MACPORTS_PREFIX}/include ${CMAKE_INCLUDE_PATH})
        MESSAGE(STATUS "Detected MacPorts installation: ${MACPORTS_PREFIX}")
    ELSE()
        UNSET(MACPORTS_PREFIX)
    ENDIF()

    IF (${DETECT_HOMEBREW} EQUAL 0)
        SET(CMAKE_LIBRARY_PATH ${HOMEBREW_PREFIX}/lib ${CMAKE_LIBRARY_PATH})
        SET(CMAKE_INCLUDE_PATH ${HOMEBREW_PREFIX}/include ${CMAKE_INCLUDE_PATH})
        MESSAGE(STATUS "Detected Homebrew installation: ${HOMEBREW_PREFIX}")
    ELSE()
        UNSET(HOMEBREW_PREFIX)
    ENDIF()

    UNSET(DETECT_HOMEBREW)
    UNSET(DETECT_MACPORTS)
ENDIF()

# Attempt to retrieve git branch and SHA1 hash of current changeset.
INCLUDE(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Set library, binary, include, share and doc paths.
SET(LIB_DIR "lib")
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND NOT APPLE )
    SET(LIB_DIR "lib64")
ENDIF()
SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")

SET(NEKTAR_BIN_DIR bin)
SET(NEKTAR_LIB_DIR ${LIB_DIR})
SET(NEKTAR_INCLUDE_DIR include/nektar++)
SET(NEKTAR_SHARE_DIR share/nektar++)
SET(NEKTAR_DOC_DIR share/doc/nektar++)
SET(NEKTAR_CMAKE_DIR ${LIB_DIR}/nektar++/cmake)

# Set ThirdParty locations.
SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)
SET(TPINC   ${CMAKE_INSTALL_PREFIX}/include/ThirdParty)
SET(TPLIB   ${CMAKE_INSTALL_PREFIX}/${NEKTAR_LIB_DIR}/nektar++/thirdparty)

# Create ThirdParty source directory if it doesn't exist already.
IF (NOT EXISTS ${TPSRC})
    FILE(MAKE_DIRECTORY ${TPSRC})
ENDIF ()

# Build shared libraries.
SET(NEKTAR_LIBRARY_TYPE "SHARED")

# Set up RPATH
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_PREFIX}/${LIB_DIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_DIR}" "${TPLIB}/")
ELSE()
    SET(CMAKE_INSTALL_RPATH "${TPLIB}")
ENDIF()

# Enable the use of @rpath in macOS install names so that we can use multiple
# third-party directories.
IF(NOT DEFINED CMAKE_MACOSX_RPATH)
    SET(CMAKE_MACOSX_RPATH 1)
ENDIF()

# Components of the Nektar++ package to build
OPTION(NEKTAR_BUILD_LIBRARY    "Build main Nektar++ libraries." ON)
OPTION(NEKTAR_BUILD_DEMOS      "Build demonstration codes."     ON)
OPTION(NEKTAR_BUILD_SOLVERS    "Build example solvers."         ON)
OPTION(NEKTAR_BUILD_UTILITIES  "Build utilities."               ON)
OPTION(NEKTAR_BUILD_UNIT_TESTS "Build unit tests."              ON)
OPTION(NEKTAR_BUILD_TESTS      "Build regression tests."        ON)
OPTION(NEKTAR_BUILD_TIMINGS    "Build benchmark timing codes."  OFF)
OPTION(NEKTAR_BUILD_PYTHON     "Build Nektar++ Python bindings" OFF)

OPTION(NEKTAR_TEST_ALL "Include full set of regression tests to this build." OFF)
OPTION(NEKTAR_TEST_USE_HOSTFILE "Use a hostfile to explicitly specify number of
slots." OFF)

OPTION(NEKTAR_ERROR_ON_WARNINGS "Use -Werror to make warnings stop compilation." ON)
MARK_AS_ADVANCED(NEKTAR_ERROR_ON_WARNINGS)

# Meshing utilities and library
IF (NOT WIN32)
    OPTION(NEKTAR_USE_MESHGEN "Build mesh generation utilities." OFF)
ENDIF()

# Build options
OPTION(NEKTAR_FULL_DEBUG "Enable Full Debugging." OFF)
MARK_AS_ADVANCED(NEKTAR_FULL_DEBUG)

IF (${CMAKE_COMPILER_IS_GNUCXX})
    OPTION(NEKTAR_ENABLE_PROFILE "Uses -pg compiler flag" OFF)
    MARK_AS_ADVANCED(NEKTAR_ENABLE_PROFILE)
ENDIF (${CMAKE_COMPILER_IS_GNUCXX})

# Memory pools
OPTION(NEKTAR_USE_MEMORY_POOLS
    "Use memory pools to accelerate memory allocation." ON)
MARK_AS_ADVANCED(NEKTAR_USE_MEMORY_POOLS)

# Thread safety
OPTION(NEKTAR_USE_THREAD_SAFETY
    "Guarantee thread safety in certain core Nektar++ classes." OFF)
MARK_AS_ADVANCED(NEKTAR_USE_THREAD_SAFETY)
IF (NEKTAR_USE_THREAD_SAFETY)
    ADD_DEFINITIONS(-DNEKTAR_USE_THREAD_SAFETY)
ENDIF()

# Disable backups
OPTION(NEKTAR_DISABLE_BACKUPS "Disable automatic backup file creation" OFF)
IF (NEKTAR_DISABLE_BACKUPS)
    ADD_DEFINITIONS(-DNEKTAR_DISABLE_BACKUPS)
ENDIF()

IF (MSVC)
    # Needed for M_PI to be visible in visual studio.
    ADD_DEFINITIONS(-D_USE_MATH_DEFINES)

    # Removes the warnings about unsafe methods such as strcpy, std::copy,
    # memcmp, etc.
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE)
ENDIF()

# Include Nektar++ common macros
INCLUDE (NektarCommon)

# Set various ThirdParty locations
OPTION(THIRDPARTY_USE_SSL "Use secure HTTP connection to download third-party files." OFF)
IF (THIRDPARTY_USE_SSL)
    SET(TPURL https://www.nektar.info/thirdparty)
ELSE()
    SET(TPURL http://www.nektar.info/thirdparty)
ENDIF()

# Find ThirdParty libraries and headers.
# --------------------------------------

# Need to find MPI first, as other packages depend on the compiler wrappers
INCLUDE (ThirdPartyMPI)
INCLUDE (ThirdPartyTinyxml)
INCLUDE (ThirdPartyMetis)
INCLUDE (ThirdPartyHDF5)
INCLUDE (ThirdPartyScotch)
INCLUDE (ThirdPartyZlib)
INCLUDE (ThirdPartyBoost)
INCLUDE (ThirdPartyPython)
INCLUDE (ThirdPartyFFTW)
INCLUDE (ThirdPartyArpack)
INCLUDE (ThirdPartyVTK)
INCLUDE (ThirdPartyOCE)
INCLUDE (ThirdPartyTriangle)
INCLUDE (ThirdPartyTetGen)
INCLUDE (ThirdPartyCCM)
INCLUDE (ThirdPartyBlasLapack)
INCLUDE (ThirdPartyCwipi)
INCLUDE (FindCFI)

INCLUDE (Doxygen)

IF( NEKTAR_USE_TINYXML_STL )
    ADD_DEFINITIONS( -DTIXML_USE_STL)
ENDIF( NEKTAR_USE_TINYXML_STL )

IF( NEKTAR_USE_MEMORY_POOLS )
    ADD_DEFINITIONS(-DNEKTAR_MEMORY_POOL_ENABLED)
ELSE( NEKTAR_USE_MEMORY_POOLS )
    REMOVE_DEFINITIONS(-DNEKTAR_MEMORY_POOL_ENABLED)
ENDIF( NEKTAR_USE_MEMORY_POOLS )

INCLUDE (ThirdPartyPETSc)

SET(Boost_USE_STATIC_LIBS OFF)
IF( WIN32 )
    # The auto-linking feature has problems with USE_STATIC_LIBS off, so we use
    # BOOST_ALL_NO_LIB to turn it off.
    # Several boost libraries headers aren't configured correctly if
    # USE_STATIC_LIBS is off, so we explicitly say they are dynamic with the
    # remaining definitions.
    ADD_DEFINITIONS(-DBOOST_ALL_NO_LIB -DBOOST_PROGRAM_OPTIONS_DYN_LINK -DBOOST_IOSTREAMS_DYN_LINK -DBOOST_THREAD_DYN_LINK)
ENDIF( )

# Set up include directories. Include work-around for -isystem for CMake (see
# bug 0010837) to avoid unnecessary compiler warnings for third-party libs.
IF (APPLE)
  SET(CMAKE_INCLUDE_SYSTEM_FLAG_C   "-isystem ")
  SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
ENDIF (APPLE)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# Build active components. Add utilities and solvers directories first, because
# that allows us to detect library dependencies automatically.
IF (NEKTAR_BUILD_LIBRARY)
    INCLUDE_DIRECTORIES(library)
ENDIF()

INCLUDE_DIRECTORIES(utilities)
ADD_SUBDIRECTORY(utilities)

INCLUDE_DIRECTORIES(solvers)
ADD_SUBDIRECTORY(solvers)

IF (NEKTAR_BUILD_LIBRARY)
    # List of Nektar++ libraries will be rebuilt every configuration.
    SET(NEKTAR++_LIBRARIES "" CACHE INTERNAL "")

    ADD_SUBDIRECTORY(library)
    INSTALL(EXPORT Nektar++Libraries DESTINATION ${NEKTAR_CMAKE_DIR} COMPONENT dev)
ENDIF (NEKTAR_BUILD_LIBRARY)

IF (NEKTAR_BUILD_TESTS)
    INCLUDE_DIRECTORIES(tests)
    ADD_SUBDIRECTORY(tests)
    IF (NEKTAR_USE_MPI)
        OPTION(NEKTAR_TEST_FORCEMPIEXEC  "Force all tests to be run through the MPI job launcher."  OFF)
    ENDIF (NEKTAR_USE_MPI)
ENDIF (NEKTAR_BUILD_TESTS)

# Compile list of definitions for Nektar++Config.cmake input file. We exclude
# vtk definitions which should not be required for external builds.
SET(NEKTAR_DEFINITIONS "")
GET_DIRECTORY_PROPERTY(
    NEKTAR_DEFINITIONS_LIST DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
FOREACH(def ${NEKTAR_DEFINITIONS_LIST})
    IF (NOT def MATCHES "^vtk")
        SET(NEKTAR_DEFINITIONS "${NEKTAR_DEFINITIONS} -D${def}")
    ENDIF()
ENDFOREACH()

# Write out Nektar++ build configuration
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/Nektar++Config.cmake.in
               ${CMAKE_BINARY_DIR}/Nektar++Config.cmake @ONLY)

# Install Nektar++ CMake configuration file
INSTALL(FILES ${CMAKE_BINARY_DIR}/Nektar++Config.cmake
    DESTINATION ${NEKTAR_CMAKE_DIR}
    COMPONENT dev)

INSTALL(CODE "FILE(GLOB tplibs ${TPDIST}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
    IF (NOT tplibs STREQUAL \"\")
        FILE(INSTALL \${tplibs} DESTINATION ${TPLIB})
    ENDIF()
")

# Install ThirdParty headers to subdirectory of ${NEKTAR_INCLUDE_DIR}
INSTALL(DIRECTORY ${TPDIST}/include/
    DESTINATION ${NEKTAR_INCLUDE_DIR}/ThirdParty
    COMPONENT dev
    OPTIONAL
    )

ADD_SUBDIRECTORY(docs)
