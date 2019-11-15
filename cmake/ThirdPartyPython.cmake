########################################################################
#
# ThirdParty configuration for Nektar++
#
# Python interfaces
#
########################################################################

IF (NEKTAR_BUILD_PYTHON)
    CMAKE_DEPENDENT_OPTION(NEKTAR_USE_PYTHON3
        "If true, prefer to use Python 3." OFF "NEKTAR_BUILD_PYTHON" OFF)

    # Set the Python3 status flag as the opposite of the Python3 flag by
    # default on first run to ensure Python is searched for.
    IF (NOT DEFINED NEKTAR_PYTHON3_STATUS)
        SET(NEKTAR_PYTHON3_STATUS NOT ${NEKTAR_USE_PYTHON3} CACHE INTERNAL "")
    ENDIF()

    IF (NOT NEKTAR_PYTHON3_STATUS STREQUAL NEKTAR_USE_PYTHON3)
        unset(PYTHON_EXECUTABLE CACHE)
        unset(PYTHON_INCLUDE_DIR CACHE)
        unset(PYTHON_LIBRARY CACHE)
        unset(PYTHON_LIBRARY_DEBUG CACHE)
        unset(BOOST_PYTHON_LIB CACHE)
        unset(BOOST_NUMPY_LIB CACHE)
        SET(NEKTAR_PYTHON3_STATUS ${NEKTAR_USE_PYTHON3} CACHE INTERNAL "")
    ENDIF()

    SET(PYTHONVER 2.7)
    IF (NEKTAR_USE_PYTHON3)
        SET(PYTHONVER 3.0)
    ENDIF()

    # Find Python
    FIND_PACKAGE(PythonInterp  ${PYTHONVER} REQUIRED)
    FIND_PACKAGE(PythonLibsNew ${PYTHONVER} REQUIRED)

    # Include headers from root directory for config file.

    # Now try to find Boost::Python. For now we are relying entirely on
    # distributed versions of this (versus trying to compile via ThirdParty)
    # because they come with various names and FindBoost is not really geared up
    # for this at present. Therefore this is done separately to avoid lots of
    # warnings and extraneous output.
    #
    # We need to try a few variants, depending on if we're doing Python 2 or
    # Python 3. Irritatingly this is all very much distriution dependent so we
    # just take our best guess at filenames. Seemingly from Boost 1.67 onwards,
    # names are just `python27` and `numpy32` but for now we have to deal with
    # this ourselves.
    STRING(REPLACE "." ";" BOOST_PYTHON_VERSION ${PYTHONLIBS_VERSION_STRING})
    LIST(GET BOOST_PYTHON_VERSION 0 BOOST_PYTHON_VERSION_MAJOR)
    LIST(GET BOOST_PYTHON_VERSION 1 BOOST_PYTHON_VERSION_MINOR)
    SET(TMP_BOOST_LIST ${NEEDED_BOOST_LIBS})

    # If we're using multi-threaded, the existing library likely has a '-mt'
    # suffix so we need to append this too.
    IF (Boost_SYSTEM_LIBRARY MATCHES "-mt")
        SET(BOOST_LIB_SUFFIX "-mt")
    ENDIF()

    # Try to find Boost::Python
    FIND_LIBRARY(BOOST_PYTHON_LIB
        NAMES boost_python-py${BOOST_PYTHON_VERSION_MAJOR}${BOOST_PYTHON_VERSION_MINOR}${BOOST_LIB_SUFFIX}
        boost_python${BOOST_PYTHON_VERSION_MAJOR}${BOOST_PYTHON_VERSION_MINOR}${BOOST_LIB_SUFFIX}
        boost_python-py${BOOST_PYTHON_VERSION_MAJOR}${BOOST_LIB_SUFFIX}
        boost_python${BOOST_PYTHON_VERSION_MAJOR}${BOOST_LIB_SUFFIX}
        boost_python${BOOST_LIB_SUFFIX}
        PATHS ${Boost_LIBRARY_DIRS})

    IF (NOT BOOST_PYTHON_LIB)
        MESSAGE(FATAL_ERROR "Could not find Boost Python installation: check it is installed in your distribution.")
    ENDIF()

    # Try to find Boost.NumPy
    FIND_LIBRARY(BOOST_NUMPY_LIB
        NAMES boost_numpy-py${BOOST_PYTHON_VERSION_MAJOR}${BOOST_PYTHON_VERSION_MINOR}${BOOST_LIB_SUFFIX}
        boost_numpy${BOOST_PYTHON_VERSION_MAJOR}${BOOST_PYTHON_VERSION_MINOR}${BOOST_LIB_SUFFIX}
        boost_numpy-py${BOOST_PYTHON_VERSION_MAJOR}${BOOST_LIB_SUFFIX}
        boost_numpy${BOOST_PYTHON_VERSION_MAJOR}${BOOST_LIB_SUFFIX}
        boost_numpy${BOOST_LIB_SUFFIX}
        PATHS ${Boost_LIBRARY_DIRS})

    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
    ADD_DEFINITIONS(-DWITH_PYTHON)

    MESSAGE(STATUS "Found Python: ${PYTHON_EXECUTABLE}")
    MESSAGE(STATUS "Found Boost.Python: ${BOOST_PYTHON_LIB}")

    # If we can't find it, pull it from git and compile it
    IF (NOT BOOST_NUMPY_LIB)
        INCLUDE(ExternalProject)
        IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            SET(WARNING_FLAGS "-Wno-cpp")
        ENDIF()
        IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            SET(WARNING_FLAGS "-Wno-#warnings")
        ENDIF()

        EXTERNALPROJECT_ADD(
            boost-numpy
            PREFIX ${TPSRC}
            URL ${TPURL}/boost-numpy_1.0.2.tar.bz2
            URL_MD5 250a517556e67f65c8837c73f419f773
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/boost-numpy
            BINARY_DIR ${TPBUILD}/boost-numpy
            TMP_DIR ${TPBUILD}/boost-numpy-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR} -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                -DCMAKE_CXX_FLAGS=${WARNING_FLAGS}
                -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
                -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DLIBRARY_TYPE=STATIC
                ${TPSRC}/boost-numpy
            )
        IF (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
            SET(BOOST_NUMPY_LIB ${TPDIST}/lib64/${CMAKE_STATIC_LIBRARY_PREFIX}boost_numpy${CMAKE_STATIC_LIBRARY_SUFFIX})
        ELSE()
            SET(BOOST_NUMPY_LIB ${TPDIST}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}boost_numpy${CMAKE_STATIC_LIBRARY_SUFFIX})
        ENDIF()
        
        INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include)
        MESSAGE(STATUS "Build Boost.NumPy: ${BOOST_NUMPY_LIB}")
    ELSE()
        MESSAGE(STATUS "Found Boost.NumPy: ${BOOST_NUMPY_LIB}")
        ADD_CUSTOM_TARGET(boost-numpy ALL)
        ADD_DEFINITIONS(-DBOOST_HAS_NUMPY)
    ENDIF()
    
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/python/setup.py.in ${CMAKE_BINARY_DIR}/setup.py)

    ADD_CUSTOM_TARGET(nekpy-install-user
        DEPENDS _MultiRegions
        COMMAND ${PYTHON_EXECUTABLE} setup.py install --user
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    ADD_CUSTOM_TARGET(nekpy-install-system
        DEPENDS _MultiRegions
        COMMAND ${PYTHON_EXECUTABLE} setup.py install
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/NekPy)
    FILE(WRITE ${CMAKE_BINARY_DIR}/NekPy/__init__.py "# Adjust dlopen flags to avoid OpenMPI issues
try:
    import DLFCN as dl
    import sys
    sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)
except ImportError:
    pass")

    MARK_AS_ADVANCED(BOOST_PYTHON_LIB)
    MARK_AS_ADVANCED(BOOST_NUMPY_LIB)
ENDIF()
