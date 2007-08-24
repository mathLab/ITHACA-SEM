SET(StdProject1DSource  Project1D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject1D StdProject1DSource)
    TARGET_LINK_LIBRARIES(StdProject1D
    optimized StdRegions debug StdRegions-g    
    )
    SET_LAPACK_LINK_LIBRARIES(StdProject1D)

SET(StdProject_Diff1DSource  Project_Diff1D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject_Diff1D StdProject_Diff1DSource)
    TARGET_LINK_LIBRARIES(StdProject_Diff1D
    optimized StdRegions debug StdRegions-g    
    )
    SET_LAPACK_LINK_LIBRARIES(StdProject_Diff1D)

SET(StdProject2DSource  Project2D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject2D StdProject2DSource)
    TARGET_LINK_LIBRARIES(StdProject2D
    optimized StdRegions debug StdRegions-g    
    )
    SET_LAPACK_LINK_LIBRARIES(StdProject2D)


