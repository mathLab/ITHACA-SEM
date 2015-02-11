MACRO(ADD_UTILITIES_EXECUTABLE name component source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})

    SET_COMMON_PROPERTIES(${name})

    TARGET_LINK_LIBRARIES(${name} MultiRegions)
    
    IF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
        SET_TARGET_PROPERTIES(${name} 
            PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup -Wl,-rpath,${CMAKE_INSTALL_PREFIX}/${LIB_DIR} -Wl,-rpath,${Boost_LIBRARY_DIRS}")
    ENDIF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
	SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${component})
    INSTALL(TARGETS ${name} 
		RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component})

ENDMACRO(ADD_UTILITIES_EXECUTABLE name component source)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library ${CMAKE_SOURCE_DIR}/utilities)

SUBDIRS(PreProcessing PostProcessing)

SET(CPACK_COMPONENTS_ALL util)

SET(CPACK_COMPONENT_UTIL_NAME "nektar++-util")
SET(CPACK_COMPONENT_UTIL_DESCRIPTION "Utilities for Nektar++")
#SET(CPACK_COMPONENT_SOLVERS_DEPENDS lib)

INCLUDE(CPack)
