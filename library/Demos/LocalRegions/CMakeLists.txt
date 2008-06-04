SET(LinkLibraries 
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g 
    optimized LibUtilities debug LibUtilities-g
)

SET(LocProject1DSource  Project1D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject1D LocProject1DSource)
TARGET_LINK_LIBRARIES(LocProject1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject1D)


SET(LocProject2DSource  Project2D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject2D LocProject2DSource)
TARGET_LINK_LIBRARIES(LocProject2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject2D)


SET(LocProject_Diff1DSource  Project_Diff1D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject_Diff1D LocProject_Diff1DSource)
TARGET_LINK_LIBRARIES(LocProject_Diff1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject_Diff1D)


