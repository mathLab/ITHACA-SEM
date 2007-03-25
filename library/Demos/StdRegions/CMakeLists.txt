SET(Project1DSource  Project1D.cpp)
ADD_NEKTAR_EXECUTABLE(Project1D Project1DSource)
    TARGET_LINK_LIBRARIES(Project1D
    optimized StdRegions debug StdRegions-g    
    )

SET(Project_Diff1DSource  Project_Diff1D.cpp)
ADD_NEKTAR_EXECUTABLE(Project_Diff1D Project_Diff1DSource)
    TARGET_LINK_LIBRARIES(Project_Diff1D
    optimized StdRegions debug StdRegions-g    
    )

IF( USE_FRAMEWORK )
    TARGET_LINK_LIBRARIES(Project1D ${FRAMEWORK_LINK_FLAGS})
    TARGET_LINK_LIBRARIES(Project_Diff1D ${FRAMEWORK_LINK_FLAGS})
ENDIF ( USE_FRAMEWORK )

