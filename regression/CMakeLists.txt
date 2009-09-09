CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(REGRESSION)

MACRO(ADD_REGRESSION_EXECUTABLE name source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})

ENDMACRO(ADD_REGRESSION_EXECUTABLE name source)

SUBDIRS(Demo Solver Auxiliary)

# Add Definitions to work Path

IF( WIN32 )
	ADD_DEFINITIONS(-DUSR_WIN_32="1")
ENDIF(WIN32)

ADD_DEFINITIONS(-DREG_PATH="${CMAKE_SOURCE_DIR}")


# CREATE MAIN REGRESSION TEST

SET(regTestSource
	./Auxiliary/regBase.cpp
	./Demo/regDemo.cpp
	./Solver/regSolver.cpp
	./main.cpp
)       

ADD_REGRESSION_EXECUTABLE(regTest
			 ${regTestSource})


# CREATE MAKE REGRESSION TEST

SET(makeRegSource
	./Auxiliary/regBase.cpp
	./Demo/regDemo.cpp
	./Solver/regSolver.cpp
	./makeReg.cpp
)

ADD_REGRESSION_EXECUTABLE(makeReg
			 ${makeRegSource})