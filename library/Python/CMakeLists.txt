CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

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


# Default install location
IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/dist CACHE PATH "" FORCE)
ENDIF()

# Find Nektar++
FIND_PACKAGE(Nektar++ REQUIRED)
ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})
INCLUDE_DIRECTORIES(SYSTEM ${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS} ${NEKTAR++_TP_LIBRARY_DIRS})
MESSAGE(STATUS "Found Nektar++: version ${NEKTAR++_VERSION}")

# Set up RPATH so that we grab Nektar++ libraries
SET(CMAKE_MACOSX_RPATH ON)
SET(CMAKE_INSTALL_RPATH ${NEKTAR++_LIBRARY_DIRS})

# Find Boost.Python
FIND_PACKAGE(Boost REQUIRED COMPONENTS python system)
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})
SET(BOOST_TMP_LIBS "${Boost_LIBRARIES}")

# Find Python
FIND_PACKAGE(PythonInterp  REQUIRED)
FIND_PACKAGE(PythonLibsNew REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_DIRS})

# Fix RPATH issue on OS X
SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF (NOT EXISTS ${TPSRC})
    FILE(MAKE_DIRECTORY ${TPSRC})
ENDIF ()

# Include headers from root directory for config file.
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# Try to find Boost.NumPy
FIND_LIBRARY(BOOST_NUMPY_LIB boost_numpy PATHS ${Boost_LIBRARY_DIRS})

# If we can't find it, pull it from git and compile it
IF (NOT BOOST_NUMPY_LIB)
    INCLUDE(ExternalProject)
    EXTERNALPROJECT_ADD(
        boost-numpy
        PREFIX ${TPSRC}
        GIT_REPOSITORY http://github.com/ndarray/Boost.NumPy.git
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/boost-numpy
        BINARY_DIR ${TPBUILD}/boost-numpy
        TMP_DIR ${TPBUILD}/boost-numpy-tmp
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR} -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            ${TPSRC}/boost-numpy
        )

    message(STATUS ${TPDIST})
    SET(BOOST_NUMPY_LIB ${TPDIST}/lib64/${CMAKE_SHARED_LIBRARY_PREFIX}boost_numpy${CMAKE_SHARED_LIBRARY_SUFFIX})
    INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include)

    # Add this directory to the RPATH.
    SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} ${TPDIST}/lib64)
ELSE()
    ADD_CUSTOM_TARGET(boost-numpy ALL)
    ADD_DEFINITIONS(-DBOOST_HAS_NUMPY)
ENDIF()

# Define a macro which will create our libraries.
MACRO(ADD_NEKPY_LIBRARY name sources)
    ADD_LIBRARY(_${name} SHARED ${${sources}})

    # Python requires a .so extension, even on OS X.
    SET_TARGET_PROPERTIES(_${name} PROPERTIES PREFIX "")
    SET_TARGET_PROPERTIES(_${name} PROPERTIES SUFFIX ".so")

    ADD_DEPENDENCIES(_${name} boost-numpy)

    # Add target link libraries.
    TARGET_LINK_LIBRARIES(_${name}
        ${BOOST_TMP_LIBS}
        ${BOOST_NUMPY_LIB}
        ${PYTHON_LIBRARIES}
        ${NEKTAR++_LIBRARIES}
        ${NEKTAR++_TP_LIBRARIES}
        )

    INSTALL(TARGETS _${name} DESTINATION
        ${CMAKE_INSTALL_PREFIX}/NekPy/${name})
ENDMACRO()

# Add libraries.
ADD_SUBDIRECTORY(LibUtilities)
ADD_SUBDIRECTORY(StdRegions)
ADD_SUBDIRECTORY(SpatialDomains)
ADD_SUBDIRECTORY(LocalRegions)
ADD_SUBDIRECTORY(MultiRegions)

# Copy lib structure.
INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/lib/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/NekPy)

# Copy examples
INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/example/
    DESTINATION ${CMAKE_INSTALL_PREFIX})
