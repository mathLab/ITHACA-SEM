######################################################
# This CMakeLists.txt is used to generate an executable
# that uses Nektar++.  Customization points for your
# projects are marked below.
######################################################

# Require at least CMake 3.1 for CMAKE_CXX_STANDARD.
CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

# TODO - Put the name of your project here.
PROJECT(ExampleSolver)

# Finds Nektar++. This may introduce other dependencies, such as Boost.
FIND_PACKAGE(Nektar++ REQUIRED)

# Require C++11
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add various definitions from Nektar++ CMake file.
ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})
INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})

# If your libraries are in a non-standard link location, you might need to add
# the link directory to the RPATH.
SET(CMAKE_INSTALL_RPATH "${NEKTAR++_LIBRARY_DIRS}")

# Set a default install location
IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/dist CACHE PATH "" FORCE)
ENDIF()

# Add a sample executable based on a driver and equation system.
ADD_EXECUTABLE(ExampleSolver
    ExampleSolver.cpp
    EquationSystems/Laplace.cpp
    )

# If your executable needs to link to any other libraries, add them here.
TARGET_LINK_LIBRARIES(ExampleSolver
    ${NEKTAR++_LIBRARIES}
    ${NEKTAR++_TP_LIBRARIES}
    )

INSTALL(TARGETS ExampleSolver DESTINATION ${CMAKE_BINARY_DIR}/dist)
