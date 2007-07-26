SET(Project1DSource  Project1D.cpp)
ADD_NEKTAR_EXECUTABLE(Project1D Project1DSource)
    TARGET_LINK_LIBRARIES(Project1D
    optimized StdRegions debug StdRegions-g    
    )
    SET_LAPACK_LINK_LIBRARIES(Project1D)

SET(Project_Diff1DSource  Project_Diff1D.cpp)
ADD_NEKTAR_EXECUTABLE(Project_Diff1D Project_Diff1DSource)
    TARGET_LINK_LIBRARIES(Project_Diff1D
    optimized StdRegions debug StdRegions-g    
    )
    SET_LAPACK_LINK_LIBRARIES(Project_Diff1D)

SET(Project2DSource  Project2D.cpp)
ADD_NEKTAR_EXECUTABLE(Project2D Project2DSource)
    TARGET_LINK_LIBRARIES(Project2D
    optimized StdRegions debug StdRegions-g    
    )
    SET_LAPACK_LINK_LIBRARIES(Project2D)


