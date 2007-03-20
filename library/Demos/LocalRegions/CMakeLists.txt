SET(Project1DSource  Project1D.cpp)
ADD_NEKTAR_EXECUTABLE(Project1D Project1DSource)
    TARGET_LINK_LIBRARIES(Project1D
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g    
    )

IF( USE_FRAMEWORK )
    TARGET_LINK_LIBRARIES(Project1D ${FRAMEWORK_LINK_FLAGS})
ENDIF ( USE_FRAMEWORK )

