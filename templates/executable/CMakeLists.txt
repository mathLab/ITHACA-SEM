######################################################
# This CMakeLists.txt is used to generate an executable
# that uses Nektar++.  Customization points for your
# projects are marked below.
######################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# TODO - Put the name of your project here.
SET(ProjectName ExampleSolver)
PROJECT(${ProjectName})

# Setup where CMake will look for modules.  You probably
# won't need to modify this.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})

# Finds Nektar++.  This may introduce other dependencies, such
# as Boost and Loki.
FIND_PACKAGE(Nektar++ REQUIRED)
ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})
INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS} ${NEKTAR++_TP_LIBRARY_DIRS})

# TODO - Add an entry here for every source file in your project. Must not be empty!
SET(SourceFileNames
    ExampleSolver.cpp
    EquationSystems/Laplace.cpp
)

# TODO - Add an entry here for every header file in your project. Must not be empty!
SET(HeaderFileNames
    EquationSystems/Laplace.h
)

ADD_EXECUTABLE(${ProjectName} ${SourceFileNames} ${HeaderFileNames})

# TODO - If your executable needs to link to any other libraries,
# add them here.
TARGET_LINK_LIBRARIES(${ProjectName}
    ${NEKTAR++_LIBRARIES}
    ${NEKTAR++_TP_LIBRARIES}
    )
