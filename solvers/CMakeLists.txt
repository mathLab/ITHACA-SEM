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

IF( NOT(CMAKE_INSTALL_PREFIX_SET) )
	SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/dist CACHE FILEPATH "Install Path" FORCE)
ENDIF()
SET(CMAKE_INSTALL_PREFIX_SET ON CACHE INTERNAL "")

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

    INCLUDE_DIRECTORIES(${SOLVERS_SOURCE_DIR})        
ENDIF( NEKTAR++_FOUND )

MACRO(ADD_SOLVER_EXECUTABLE name source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})

    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)

    IF( NEKTAR++_FOUND )
        # TODO - If your executable needs to link to any other libraries,  
        # add them here.
        TARGET_LINK_LIBRARIES(${name}
	    ${NEKTAR++_LIBRARIES}
	    ${Boost_THREAD_LIBRARY} 
            ${Boost_IOSTREAMS_LIBRARY} 
	    ${Boost_ZLIB_LIBRARY}  
            ${Boost_DATE_TIME_LIBRARY}      
            optimized tinyxml debug tinyxml)

        SET_LAPACK_LINK_LIBRARIES(${name})

    ENDIF( NEKTAR++_FOUND )

    INSTALL(TARGETS ${name} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

ENDMACRO(ADD_SOLVER_EXECUTABLE name source)

# Generate list of available solvers
FILE(GLOB dir_list "*")
LIST(REMOVE_ITEM dir_list "${CMAKE_SOURCE_DIR}/Auxiliary")
LIST(REMOVE_ITEM dir_list "${CMAKE_SOURCE_DIR}/.svn")
LIST(REMOVE_ITEM dir_list "${CMAKE_SOURCE_DIR}/builds")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)    
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})
