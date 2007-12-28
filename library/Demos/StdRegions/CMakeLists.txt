SET(StdProject1DSource  Project1D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject1D StdProject1DSource)
    TARGET_LINK_LIBRARIES(StdProject1D
    optimized StdRegions debug StdRegions-g    
    )
#     SET_LAPACK_LINK_LIBRARIES(StdProject1D)

SET(StdProject_Diff1DSource  Project_Diff1D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject_Diff1D StdProject_Diff1DSource)
    TARGET_LINK_LIBRARIES(StdProject_Diff1D
    optimized StdRegions debug StdRegions-g    
    )
#     SET_LAPACK_LINK_LIBRARIES(StdProject_Diff1D)

SET(StdProject2DSource  Project2D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject2D StdProject2DSource)
    TARGET_LINK_LIBRARIES(StdProject2D
    optimized StdRegions debug StdRegions-g    
    )
#     SET_LAPACK_LINK_LIBRARIES(StdProject2D)

SET(StdTetExpDemoSource StdTetExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(StdTetExpDemo StdTetExpDemoSource)
    TARGET_LINK_LIBRARIES(StdTetExpDemo
    optimized StdRegions debug StdRegions-g    
    )

SET(StdHexExpDemoSource StdHexExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(StdHexExpDemo StdHexExpDemoSource)
    TARGET_LINK_LIBRARIES(StdHexExpDemo
    optimized StdRegions debug StdRegions-g    
    )

SET(StdTetExp_Diff3DSource StdTetExp_Diff3D.cpp)
ADD_NEKTAR_EXECUTABLE(StdTetExp_Diff3D StdTetExp_Diff3DSource)
    TARGET_LINK_LIBRARIES(StdTetExp_Diff3D
    optimized StdRegions debug StdRegions-g    
    )    
    
SET(StdPrismExpDemoSource StdPrismExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(StdPrismExpDemo StdPrismExpDemoSource)
    TARGET_LINK_LIBRARIES(StdPrismExpDemo
    optimized StdRegions debug StdRegions-g    
    )
    
SET(StdPrismExp_Diff3DSource StdPrismExp_Diff3D.cpp)
ADD_NEKTAR_EXECUTABLE(StdPrismExp_Diff3D StdPrismExp_Diff3DSource)
    TARGET_LINK_LIBRARIES(StdPrismExp_Diff3D
    optimized StdRegions debug StdRegions-g    
    ) 

    

# SET(StdProject_Diff2DSource  Project_Diff2D.cpp)
# ADD_NEKTAR_EXECUTABLE(StdProject_Diff2D StdProject_Diff2DSource)
#     TARGET_LINK_LIBRARIES(StdProject_Diff2D
#     optimized StdRegions debug StdRegions-g    
#     )
