######################################################
# This CMakeLists.txt is used to generate an utilities
# project.  Customization points for your
# projects are marked below.
######################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(UTILITIES)

IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
    SET(LIB_DIR "lib64")
ENDIF()

IF( NOT(CMAKE_INSTALL_PREFIX_SET) )
    SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/dist CACHE FILEPATH "Install   Path" FORCE)
ENDIF()
SET(CMAKE_INSTALL_PREFIX_SET ON CACHE INTERNAL "")

# Setup where CMake will look for modules.  You probably 
# won't need to modify this.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})

# Set up RPATH
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
                                "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
       SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Look for Nektar++ in the standard location (Nektar++/builds) 
# Users can modify this directly in CMake if they used a non-standard location.
SET(Nektar++_DIR "${CMAKE_SOURCE_DIR}/../builds/dist" CACHE PATH "")

# Finds Nektar++.  This may introduce other dependencies, such 
# as Boost and Loki.
FIND_PACKAGE(Nektar++ REQUIRED)

IF( NEKTAR++_FOUND )
    INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS})
    LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS})
    ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})

    INCLUDE_DIRECTORIES(${UTILITIES_SOURCE_DIR})

    IF(NEKTAR_USE_MPI)
        INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
    ENDIF(NEKTAR_USE_MPI)
ENDIF( NEKTAR++_FOUND )


MACRO(ADD_UTILITIES_EXECUTABLE name component source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})

    IF( NEKTAR++_FOUND )
        SET_COMMON_PROPERTIES(${name})
        
        # TODO - If your executable needs to link to any other libraries,  
        # add them here.
        TARGET_LINK_LIBRARIES(${name}
	        ${NEKTAR++_LIBRARIES}
	        ${Boost_THREAD_LIBRARY}
            ${Boost_IOSTREAMS_LIBRARY}
	        ${Boost_ZLIB_LIBRARY}
            ${Boost_DATE_TIME_LIBRARY}
            optimized tinyxml debug tinyxml
        )

        IF( ${CMAKE_SYSTEM} MATCHES "Linux.*" )
            TARGET_LINK_LIBRARIES(${name} optimized rt debug rt)
        ENDIF()
        
        SET_LAPACK_LINK_LIBRARIES(${name})

        IF (NEKTAR_USE_MPI)
            TARGET_LINK_LIBRARIES(${name} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
        ENDIF (NEKTAR_USE_MPI)
    ENDIF( NEKTAR++_FOUND )
    
    IF( ${CMAKE_SYSTEM} MATCHES "Linux.*")
        TARGET_LINK_LIBRARIES(${name} optimized rt debug rt)
    ENDIF()

    IF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
        SET_TARGET_PROPERTIES(${name} 
            PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup -Wl,-rpath,${CMAKE_INSTALL_PREFIX}/${LIB_DIR} -Wl,-rpath,${Boost_LIBRARY_DIRS}")
    ENDIF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
	
    INSTALL(TARGETS ${name} 
		RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component})

ENDMACRO(ADD_UTILITIES_EXECUTABLE name component source)


MACRO(SET_COMMON_PROPERTIES name)

    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)

ENDMACRO(SET_COMMON_PROPERTIES name)


SUBDIRS( PreProcessing  PostProcessing)

SET(CPACK_COMPONENTS_ALL util)

SET(CPACK_COMPONENT_UTIL_NAME "nektar++-util")
SET(CPACK_COMPONENT_UTIL_DESCRIPTION "Utilities for Nektar++")
#SET(CPACK_COMPONENT_SOLVERS_DEPENDS lib)

INCLUDE(CPack)
