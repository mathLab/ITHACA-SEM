MACRO(ADD_UTILITIES_EXECUTABLE name component source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})

    SET_COMMON_PROPERTIES(${name})
        
    # TODO - If your executable needs to link to any other libraries,  
    # add them here.
    TARGET_LINK_LIBRARIES(${name}
	        ${NEKTAR++_LIBRARIES}
	        ${Boost_THREAD_LIBRARY}
            ${Boost_IOSTREAMS_LIBRARY}
	        ${ZLIB_LIBRARY}
            ${Boost_DATE_TIME_LIBRARY}
            ${Boost_PROGRAM_OPTIONS_LIBRARY} 
            optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
    )
        
    SET_LAPACK_LINK_LIBRARIES(${name})
            
    IF (NEKTAR_USE_MPI)
        TARGET_LINK_LIBRARIES(${name} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    ENDIF (NEKTAR_USE_MPI)
    
    IF (NEKTAR_USE_VTK)
        TARGET_LINK_LIBRARIES(${name} vtkCommon vtkIO vtkGraphics)
    ENDIF (NEKTAR_USE_VTK)
    
    IF( ${CMAKE_SYSTEM} MATCHES "Linux.*")
        TARGET_LINK_LIBRARIES(${name} optimized rt debug rt)
    ENDIF()

    IF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
        SET_TARGET_PROPERTIES(${name} 
            PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup -Wl,-rpath,${CMAKE_INSTALL_PREFIX}/${LIB_DIR} -Wl,-rpath,${Boost_LIBRARY_DIRS}")
    ENDIF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
	SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${component})
    INSTALL(TARGETS ${name} 
		RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component})

ENDMACRO(ADD_UTILITIES_EXECUTABLE name component source)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library ${CMAKE_SOURCE_DIR}/utilities)

SUBDIRS( PreProcessing  PostProcessing)

SET(CPACK_COMPONENTS_ALL util)

SET(CPACK_COMPONENT_UTIL_NAME "nektar++-util")
SET(CPACK_COMPONENT_UTIL_DESCRIPTION "Utilities for Nektar++")
#SET(CPACK_COMPONENT_SOLVERS_DEPENDS lib)

INCLUDE(CPack)
