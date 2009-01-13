######################################################
# This CMakeLists.txt is used to generate an executable
# that uses Nektar++.  Customization points for your 
# projects are marked below.
######################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(SOLVERS)

# Setup where CMake will look for modules.  You probably 
# won't need to modify this.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})

# Finds Nektar++.  This may introduce other dependencies, such 
# as Boost and Loki.
FIND_PACKAGE(Nektar++ REQUIRED)

IF( NEKTAR++_FOUND )
    INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS})
    LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS})
    ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})

    INCLUDE_DIRECTORIES(${SOLVERS_SOURCE_DIR})        
ENDIF( NEKTAR++_FOUND )

MACRO(ADD_SOLVER_EXECUTABLE name source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})

    IF( NEKTAR++_FOUND )
        # TODO - If your executable needs to link to any other libraries,  
        # add them here.
        TARGET_LINK_LIBRARIES(${name}
	    ${NEKTAR++_LIBRARIES}
	    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
            optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
	    optimized ${Boost_ZLIB_LIBRARY_RELEASE} debug ${Boost_ZLIB_LIBRARY_DEBUG}
            optimized ${Boost_DATE_TIME_LIBRARY_RELEASE} debug ${Boost_DATE_TIME_LIBRARY_DEBUG}     
            optimized tinyxml debug tinyxml)

        SET_LAPACK_LINK_LIBRARIES(${name})

    ENDIF( NEKTAR++_FOUND )
ENDMACRO(ADD_SOLVER_EXECUTABLE name source)

MACRO(SET_COMMON_PROPERTIES name)

    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)


ENDMACRO(SET_COMMON_PROPERTIES name)

#SUBDIRS(AdvectionDiffusionReactionSolver ShallowWaterSolver EulerSolver)
SUBDIRS(AdvectionDiffusionReactionSolver)
