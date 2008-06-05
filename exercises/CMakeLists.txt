######################################################
# This CMakeLists.txt is used to generate an executable
# that uses Nektar++.  Customization points for your 
# projects are marked below.
######################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(exercises)

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
ENDIF( NEKTAR++_FOUND )

MACRO(ADD_EXERCISE_EXECUTABLE name source)
    ADD_EXECUTABLE(${name} ${source})

    IF( NEKTAR++_FOUND )
        # TODO - If your executable needs to link to any other libraries,  
        # add them here.
        TARGET_LINK_LIBRARIES(${name}
	    ${NEKTAR++_LIBRARIES}
	    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
            optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
	    optimized ${Boost_ZLIB_LIBRARY_RELEASE} debug ${Boost_ZLIB_LIBRARY_DEBUG}
            optimized tinyxml debug tinyxml)

        SET_LAPACK_LINK_LIBRARIES(${name})
    ENDIF( NEKTAR++_FOUND )
ENDMACRO(ADD_EXERCISE_EXECUTABLE name source)


SET(Tutorial1Source Tutorial1.cpp)
SET(Tutorial2Source Tutorial2.cpp)
SET(Tutorial3Source Tutorial3.cpp)
ADD_EXERCISE_EXECUTABLE(Tutorial1 ${Tutorial1Source})
ADD_EXERCISE_EXECUTABLE(Tutorial2 ${Tutorial2Source})
ADD_EXERCISE_EXECUTABLE(Tutorial3 ${Tutorial3Source})

