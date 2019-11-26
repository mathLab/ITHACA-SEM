########################################################################
#
# ThirdParty configuration for Nektar++
#
# Boost
#
########################################################################

#If the user has not set BOOST_ROOT, look in a couple common places first.
MESSAGE(STATUS "Searching for Boost:")

# Minimum version and boost libraries required
SET(MIN_VER "1.56.0")
SET(NEEDED_BOOST_LIBS thread iostreams filesystem system program_options regex)

SET(Boost_NO_BOOST_CMAKE ON)
IF( BOOST_ROOT )
    SET(Boost_NO_SYSTEM_PATHS ON)
    FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS})
ELSE ()
    SET(TEST_ENV1 $ENV{BOOST_HOME})
    SET(TEST_ENV2 $ENV{BOOST_DIR})
    IF (DEFINED TEST_ENV1)
        SET(BOOST_ROOT $ENV{BOOST_HOME})
        SET(Boost_NO_SYSTEM_PATHS ON)
        FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS} )
    ELSEIF (DEFINED TEST_ENV2)
        SET(BOOST_ROOT $ENV{BOOST_DIR})
        SET(Boost_NO_SYSTEM_PATHS ON)
        FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS} )
    ELSE ()
        SET(BOOST_ROOT ${TPDIST})
        FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS} )
    ENDIF()
ENDIF()

# Check what we found and determine if we need to build boost
FOREACH(FOUND_VAR ${NEEDED_BOOST_LIBS})
    STRING(TOUPPER ${FOUND_VAR} FOUND_VAR_UPPER)
    IF (Boost_${FOUND_VAR_UPPER}_FOUND)
        MESSAGE(STATUS "-- Found Boost ${FOUND_VAR} library: "
                "${Boost_${FOUND_VAR_UPPER}_LIBRARY}")
    ELSE ()
        MESSAGE(STATUS "-- Pre-installed Boost ${FOUND_VAR} library not found")
    ENDIF()
ENDFOREACH()

IF (NOT Boost_FOUND)
    SET(BUILD_BOOST ON)
ELSE()
    SET(BUILD_BOOST OFF)
ENDIF ()


OPTION(THIRDPARTY_BUILD_BOOST "Build Boost libraries" ${BUILD_BOOST})
SET(Boost_USE_MULTITHREADED ON CACHE BOOL
    "Search for multithreaded boost libraries")
MARK_AS_ADVANCED(Boost_USE_MULTITHREADED)


IF (WIN32)
    ADD_DEFINITIONS("-DBOOST_ALL_NO_LIB")
ENDIF()

IF (THIRDPARTY_BUILD_BOOST)
    INCLUDE(ExternalProject)

    # Only build the libraries we need
    FOREACH(boostlib ${NEEDED_BOOST_LIBS})
 	LIST(APPEND BOOST_LIB_LIST --with-${boostlib})
    ENDFOREACH()

    IF (NOT WIN32)
        # We need -fPIC for 64-bit builds
        IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
            SET(BOOST_FLAGS cxxflags=-fPIC cflags=-fPIC linkflags=-fPIC)
        ENDIF ()
    ENDIF()

    # Build Boost: first need to select toolset. Some will have specific
    # versions.
    SET(TOOLSET_VERSION "")
    IF (APPLE)
        # macOS should have the darwin toolset regardless of gcc/clang.
        SET(TOOLSET darwin)
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        SET(TOOLSET msvc)
        IF (MSVC_VERSION EQUAL 1600)
            SET(TOOLSET_VERSION 10.0) # Visual Studio 2010
        ELSEIF (MSVC_VERSION EQUAL 1700)
            SET(TOOLSET_VERSION 11.0) # Visual Studio 2012
        ELSEIF (MSVC_VERSION EQUAL 1800)
            SET(TOOLSET_VERSION 12.0) # Visual Studio 2013
        ELSEIF (MSVC_VERSION EQUAL 1900)
            SET(TOOLSET_VERSION 14.0) # Visual Studio 2015
        ELSEIF (MSVC_VERSION GREATER 1909 AND MSVC_VERSION LESS 1920)
            SET(TOOLSET_VERSION 14.1) # Visual Studio 2017
        ELSEIF (MSVC_VERSION GREATER 1919 AND MSVC_VERSION LESS 1930)
            SET(TOOLSET_VERSION 14.2) # Visual Studio 2019
        ENDIF()
    ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Cray")
        SET(TOOLSET cray)
    ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        SET(TOOLSET intel)
    ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        SET(TOOLSET gcc)
        SET(TOOLSET_VERSION ${CMAKE_CXX_COMPILER_VERSION})
    ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        SET(TOOLSET clang)
        SET(TOOLSET_VERSION ${CMAKE_CXX_COMPILER_VERSION})
    ELSE()
        MESSAGE(STATUS "Unknown compiler for boost, assuming gcc toolset")
        SET(TOOLSET gcc)
    ENDIF()

    IF (TOOLSET_VERSION STREQUAL "")
        SET(TOOLSET_CMDLINE ${TOOLSET})
    ELSE()
        SET(TOOLSET_CMDLINE ${TOOLSET}-${TOOLSET_VERSION})
    ENDIF()

    IF (NOT WIN32)
        EXTERNALPROJECT_ADD(
            boost
            PREFIX ${TPSRC}
            URL ${TPURL}/boost_1_57_0.tar.bz2
            URL_MD5 "1be49befbdd9a5ce9def2983ba3e7b76"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/boost
            BINARY_DIR ${TPBUILD}/boost
            TMP_DIR ${TPBUILD}/boost-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ./bootstrap.sh
            BUILD_COMMAND NO_BZIP2=1 ./b2
                variant=release
                link=shared
                include=${TPDIST}/include
                cxxflags="-w"
                linkflags="-L${TPDIST}/lib"
                ${BOOST_FLAGS} ${BOOST_LIB_LIST}
		--prefix=${TPDIST}
                --layout=system toolset=${TOOLSET_CMDLINE} install
            INSTALL_COMMAND ""
            )
    ELSE ()
        IF (CMAKE_CL_64)
            SET(ADDRESS_MODEL 64)
        ELSE()
            SET(ADDRESS_MODEL 32)
        ENDIF()
        EXTERNALPROJECT_ADD(
            boost
            PREFIX ${TPSRC}
            URL ${TPURL}/boost_1_57_0.tar.bz2
            URL_MD5 "1be49befbdd9a5ce9def2983ba3e7b76"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/boost
            BINARY_DIR ${TPBUILD}/boost
            TMP_DIR ${TPBUILD}/boost-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND call bootstrap.bat
            BUILD_COMMAND b2 variant=release
                toolset=${TOOLSET_CMDLINE}
                address-model=${ADDRESS_MODEL}
                link=shared
                runtime-link=shared
                -s NO_BZIP2=1
                -s ZLIB_BINARY=zlib
                -s ZLIB_INCLUDE=${TPDIST}/include
                -s ZLIB_LIBPATH=${TPDIST}/lib
                ${BOOST_LIB_LIST}
                --layout=system
                --prefix=${TPDIST} install
            INSTALL_COMMAND ""
            )
    ENDIF()

    IF (APPLE)
        EXTERNALPROJECT_ADD_STEP(boost patch-install-path
            COMMAND sed -i ".bak" "s|-install_name \"|&${TPDIST}/lib/|"
                ${TPBUILD}/boost/tools/build/src/tools/darwin.jam
            DEPENDERS build
            DEPENDEES download)
    ENDIF (APPLE)

    # Write to jamfile to use appropriate toolset.
    SET(cmd_string "using ${TOOLSET} : ${TOOLSET_VERSION}")
    SET(cmd_string "${cmd_string} : ${CMAKE_CXX_COMPILER} $<SEMICOLON>")

    IF (UNIX)
	EXTERNALPROJECT_ADD_STEP(boost conf-project-conf
            COMMAND cmake -E echo "${cmd_string}" >
                ${TPBUILD}/boost/tools/build/src/user-config.jam
            DEPENDERS build
            DEPENDEES configure)
    ENDIF()

    # If building ThirdParty zlib, force zlib build before boost
    IF (THIRDPARTY_BUILD_ZLIB)
        ADD_DEPENDENCIES(boost zlib-1.2.7)
    ENDIF(THIRDPARTY_BUILD_ZLIB)

    # Set up CMake variables
    FOREACH(BOOSTLIB ${NEEDED_BOOST_LIBS})
        STRING(TOUPPER ${BOOSTLIB} BOOSTLIB_UPPER)
        THIRDPARTY_LIBRARY(Boost_${BOOSTLIB_UPPER}_LIBRARY
            SHARED boost_${BOOSTLIB} DESCRIPTION "Boost ${BOOSTLIB} library")
        MARK_AS_ADVANCED(Boost_${BOOSTLIB_UPPER}_LIBRARY)
        LIST(APPEND Boost_LIBRARIES ${Boost_${BOOSTLIB_UPPER}_LIBRARY})
    ENDFOREACH()

    SET(Boost_INCLUDE_DIRS ${TPSRC}/dist/include)
    SET(Boost_CONFIG_INCLUDE_DIR ${TPINC})
    SET(Boost_LIBRARY_DIRS ${TPSRC}/dist/lib)
    SET(Boost_CONFIG_LIBRARY_DIR ${TPLIB})

    INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include)

    STRING(REPLACE ";" ", " NEEDED_BOOST_LIBS_STRING "${NEEDED_BOOST_LIBS}")
    MESSAGE(STATUS "Build boost libs: ${NEEDED_BOOST_LIBS_STRING}")
ELSE (THIRDPARTY_BUILD_BOOST)
    ADD_CUSTOM_TARGET(boost ALL)
    IF (BOOST_THREAD_LIBRARY)
        MARK_AS_ADVANCED(BOOST_THREAD_LIBRARY)
    ENDIF (BOOST_THREAD_LIBRARY)
    SET(Boost_CONFIG_INCLUDE_DIR ${Boost_INCLUDE_DIRS})
    SET(Boost_CONFIG_LIBRARY_DIR ${Boost_LIBRARY_DIRS})
ENDIF (THIRDPARTY_BUILD_BOOST)

INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})
