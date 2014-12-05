MACRO(ADD_SOLVER_EXECUTABLE name component source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})
    
    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)

    SET_COMMON_PROPERTIES(${name})


    TARGET_LINK_LIBRARIES(${name} SolverUtils)

    IF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
        SET_TARGET_PROPERTIES(${name} 
            PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup -Wl,-rpath,${CMAKE_INSTALL_PREFIX}/${LIB_DIR} -Wl,-rpath,${Boost_LIBRARY_DIRS}")
    ENDIF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${component})

	INSTALL(TARGETS ${name} 
		RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component} OPTIONAL)

ENDMACRO(ADD_SOLVER_EXECUTABLE name component source)


# Generate list of available subdirectories
FILE(GLOB dir_list "*")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)    
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})


SET(CPACK_COMPONENTS_ALL solvers)

SET(CPACK_COMPONENT_SOLVERS_NAME "nektar++-solvers")
SET(CPACK_COMPONENT_SOLVERS_DESCRIPTION "Solvers for Nektar++")
#SET(CPACK_COMPONENT_SOLVERS_DEPENDS lib)

INCLUDE(CPack)
