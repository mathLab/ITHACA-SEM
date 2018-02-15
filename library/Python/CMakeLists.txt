# Find Python
FIND_PACKAGE(PythonInterp  REQUIRED)
FIND_PACKAGE(PythonLibsNew REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_DIRS})

# Include headers from root directory for config file.

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
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_PYTHON_LIBRARY}
        ${BOOST_NUMPY_LIB}
        ${PYTHON_LIBRARIES}
        ${name})

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
INSTALL(DIRECTORY lib/ DESTINATION ${CMAKE_INSTALL_PREFIX}/NekPy)

# Copy examples
INSTALL(DIRECTORY example/ DESTINATION ${CMAKE_INSTALL_PREFIX})
