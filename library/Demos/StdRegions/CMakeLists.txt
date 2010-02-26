SET(LinkLibraries 
    optimized StdRegions debug StdRegions-g
    optimized LibUtilities debug LibUtilities-g
)
	  
SET(StdProject1DSource  StdProject1D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject1D StdProject1DSource)
TARGET_LINK_LIBRARIES(StdProject1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdProject1D)

SET(StdProject_Diff1DSource  StdProject_Diff1D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject_Diff1D StdProject_Diff1DSource)
TARGET_LINK_LIBRARIES(StdProject_Diff1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdProject_Diff1D)

SET(StdProject2DSource  StdProject2D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject2D StdProject2DSource)
TARGET_LINK_LIBRARIES(StdProject2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdProject2D)

SET(StdProject_Diff2DSource  StdProject_Diff2D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject_Diff2D StdProject_Diff2DSource)
TARGET_LINK_LIBRARIES(StdProject_Diff2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdProject_Diff2D)

SET(StdProject3DSource  StdProject3D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject3D StdProject3DSource)
TARGET_LINK_LIBRARIES(StdProject3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdProject3D)

SET(StdProject_Diff3DSource  StdProject_Diff3D.cpp)
ADD_NEKTAR_EXECUTABLE(StdProject_Diff3D StdProject_Diff3DSource)
TARGET_LINK_LIBRARIES(StdProject_Diff3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdProject_Diff3D)

SET(StdTetExpDemoSource StdTetExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(StdTetExpDemo StdTetExpDemoSource)
TARGET_LINK_LIBRARIES(StdTetExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdTetExpDemo)

SET(StdHexExpDemoSource StdHexExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(StdHexExpDemo StdHexExpDemoSource)
TARGET_LINK_LIBRARIES(StdHexExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdHexExpDemo)

SET(StdTetExp_Diff3DSource StdTetExp_Diff3D.cpp)
ADD_NEKTAR_EXECUTABLE(StdTetExp_Diff3D StdTetExp_Diff3DSource)
TARGET_LINK_LIBRARIES(StdTetExp_Diff3D ${LinkLibraries})  
SET_LAPACK_LINK_LIBRARIES(StdTetExp_Diff3D)   

SET(StdHexExp_Diff3DSource StdHexExp_Diff3D.cpp)
ADD_NEKTAR_EXECUTABLE(StdHexExp_Diff3D StdHexExp_Diff3DSource)
TARGET_LINK_LIBRARIES(StdHexExp_Diff3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdHexExp_Diff3D)


SET(StdPrismExpDemoSource StdPrismExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(StdPrismExpDemo StdPrismExpDemoSource)
TARGET_LINK_LIBRARIES(StdPrismExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdPrismExpDemo)
   
SET(StdPrismExp_Diff3DSource StdPrismExp_Diff3D.cpp)
ADD_NEKTAR_EXECUTABLE(StdPrismExp_Diff3D StdPrismExp_Diff3DSource)
TARGET_LINK_LIBRARIES(StdPrismExp_Diff3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdPrismExp_Diff3D)
    
SET(StdPyrExpDemoSource StdPyrExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(StdPyrExpDemo StdPyrExpDemoSource)
TARGET_LINK_LIBRARIES(StdPyrExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(StdPyrExpDemo)
    
SET(StdPyrExp_Diff3DSource StdPyrExp_Diff3D.cpp)
ADD_NEKTAR_EXECUTABLE(StdPyrExp_Diff3D StdPyrExp_Diff3DSource)
TARGET_LINK_LIBRARIES(StdPyrExp_Diff3D ${LinkLibraries})     
SET_LAPACK_LINK_LIBRARIES(StdPyrExp_Diff3D)
    
