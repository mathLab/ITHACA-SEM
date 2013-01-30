IF (NEKTAR_USE_VTK)
    SET(FIBRETONEKTAR_SOURCES FibreToNektar.cpp)
    ADD_EXECUTABLE(FibreToNektar ${FIBRETONEKTAR_SOURCES})
    SET_COMMON_PROPERTIES(FibreToNektar)
    TARGET_LINK_LIBRARIES(FibreToNektar
            ${NEKTAR++_LIBRARIES}
            ${Boost_THREAD_LIBRARY}
            ${Boost_IOSTREAMS_LIBRARY}
            ${Boost_ZLIB_LIBRARY}
            ${Boost_DATE_TIME_LIBRARY}
            ${Boost_PROGRAM_OPTIONS_LIBRARY} 
            optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
    )
    TARGET_LINK_LIBRARIES(FibreToNektar ${TINYXML_LIB})
    TARGET_LINK_LIBRARIES(FibreToNektar vtkCommon vtksys vtkViews vtkWidgets vtkRendering vtkIO)
    
    SET_LAPACK_LINK_LIBRARIES(FibreToNektar)
            
    IF (NEKTAR_USE_MPI)
        TARGET_LINK_LIBRARIES(FibreToNektar ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    ENDIF (NEKTAR_USE_MPI)
    
    INSTALL(TARGETS FibreToNektar 
		RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT solvers OPTIONAL)
ENDIF (NEKTAR_USE_VTK)
