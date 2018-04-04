########################################################################
#
# ThirdParty configuration for Nektar++
#
# Boost
#
########################################################################

# Define a macro that checks for common Boost environment variables.
MACRO(NektarFindBoost)
    SET(Boost_DEBUG 0)
    SET(Boost_NO_BOOST_CMAKE ON)
    IF (BOOST_ROOT)
        SET(Boost_NO_SYSTEM_PATHS ON)
        FIND_PACKAGE(Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS})
    ELSE ()
        SET(TEST_ENV1 $ENV{BOOST_HOME})
        SET(TEST_ENV2 $ENV{BOOST_DIR})
        IF (DEFINED TEST_ENV1)
            SET(BOOST_ROOT $ENV{BOOST_HOME})
            SET(Boost_NO_SYSTEM_PATHS ON)
            FIND_PACKAGE(Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS})
        ELSEIF (DEFINED TEST_ENV2)
            SET(BOOST_ROOT $ENV{BOOST_DIR})
            SET(Boost_NO_SYSTEM_PATHS ON)
            FIND_PACKAGE(Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS})
        ELSE ()
            SET(BOOST_ROOT ${TPDIST})
            FIND_PACKAGE(Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS})
        ENDIF()
    ENDIF()
ENDMACRO()

#If the user has not set BOOST_ROOT, look in a couple common places first.
MESSAGE(STATUS "Searching for Boost:")

# Minimum version and boost libraries required
SET(MIN_VER "1.56.0")
SET(NEEDED_BOOST_LIBS thread iostreams filesystem system program_options regex)

IF (NEKTAR_BUILD_PYTHON)
    # We need to try a few variants, depending on if we're doing Python 2 or
    # Python 3. Irritatingly this is all very much distriution dependent so we
    # just take our best guess at filenames.
    STRING(REPLACE "." ";" BOOST_PYTHON_VERSION ${PYTHONLIBS_VERSION_STRING})
    LIST(GET BOOST_PYTHON_VERSION 0 BOOST_PYTHON_VERSION_MAJOR)
    LIST(GET BOOST_PYTHON_VERSION 1 BOOST_PYTHON_VERSION_MINOR)
    SET(TMP_BOOST_LIST ${NEEDED_BOOST_LIBS})

    # Now try different versions of the string.
    SET(NAMES_TO_TRY python-py${BOOST_PYTHON_VERSION_MAJOR}${BOOST_PYTHON_VERSION_MINOR} python-py${BOOST_PYTHON_VERSION_MAJOR} python${BOOST_PYTHON_VERSION_MAJOR} python)

    FOREACH(BPLIB ${NAMES_TO_TRY})
        SET(NEEDED_BOOST_LIBS ${TMP_BOOST_LIST} ${BPLIB})
        NektarFindBoost()
        STRING(TOUPPER ${BPLIB} BPLIB_UPPER)

        IF (Boost_${BPLIB_UPPER}_FOUND)
            break()
        ENDIF()
    ENDFOREACH()
ELSE()
    NektarFindBoost()
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

    # Build Boost
    IF (APPLE)
        SET(TOOLSET darwin)
    ELSEIF (WIN32)
        IF (MSVC10)
            SET(TOOLSET msvc-10.0)
        ELSEIF (MSVC11)
            SET(TOOLSET msvc-11.0)
        ELSEIF (MSVC12)
            SET(TOOLSET msvc-12.0)
        ELSEIF (MSVC14)
            SET(TOOLSET msvc-14.0)
        ENDIF()
    ELSE()
        SET(TOOLSET gcc)
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
            CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ./bootstrap.sh --prefix=${TPDIST}
            BUILD_COMMAND NO_BZIP2=1 ./b2
                variant=release
                link=shared
                include=${TPDIST}/include
                linkflags="-L${TPDIST}/lib"
                ${BOOST_FLAGS} ${BOOST_LIB_LIST}
                --layout=system toolset=${TOOLSET} install
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
                toolset=${TOOLSET}
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
            COMMAND sed -i ".bak" "s|-install_name \"|&${TPDIST}/lib/|" ${TPBUILD}/boost/tools/build/src/tools/darwin.jam
            DEPENDERS build
            DEPENDEES download)
    ENDIF (APPLE)

    # If building ThirdParty zlib, force zlib build before boost
    IF (THIRDPARTY_BUILD_ZLIB)
        ADD_DEPENDENCIES(boost zlib-1.2.7)
    ENDIF(THIRDPARTY_BUILD_ZLIB)

    # Set up CMake variables
    FOREACH(BOOSTLIB ${NEEDED_BOOST_LIBS})
        STRING(TOUPPER ${BOOSTLIB} BOOSTLIB_UPPER)
        THIRDPARTY_LIBRARY(Boost_${BOOSTLIB_UPPER}_LIBRARY
            SHARED boost_${BOOSTLIB} DESCRIPTION "Boost ${BOOSTLIB} library")
        THIRDPARTY_LIBRARY(Boost_${BOOSTLIB_UPPER}_LIBRARY_DEBUG
            SHARED boost_${BOOSTLIB} DESCRIPTION "Boost ${BOOSTLIB} library, debug")
        THIRDPARTY_LIBRARY(Boost_${BOOSTLIB_UPPER}_LIBRARY_RELEASE
            SHARED boost_${BOOSTLIB} DESCRIPTION "Boost ${BOOSTLIB} library, release")
        MARK_AS_ADVANCED(Boost_${BOOSTLIB_UPPER}_LIBRARY)
        MARK_AS_ADVANCED(Boost_${BOOSTLIB_UPPER}_LIBRARY_DEBUG)
        MARK_AS_ADVANCED(Boost_${BOOSTLIB_UPPER}_LIBRARY_RELEASE)
        LIST(APPEND Boost_LIBRARIES ${Boost_${BOOSTLIB_UPPER}_LIBRARY})
    ENDFOREACH()

    SET(Boost_INCLUDE_DIRS ${TPSRC}/dist/include)
    SET(Boost_CONFIG_INCLUDE_DIR ${TPINC})
    SET(Boost_LIBRARY_DIRS ${TPSRC}/dist/lib)
    SET(Boost_CONFIG_LIBRARY_DIR ${TPLIB})

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

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
