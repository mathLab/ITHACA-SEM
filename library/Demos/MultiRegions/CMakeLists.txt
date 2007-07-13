SET(ProjectLoc1DSource  ProjectLoc1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc1D ProjectLoc1DSource)
    TARGET_LINK_LIBRARIES(ProjectLoc1D
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g    
    optimized tinyxml debug tinyxml
    )

SET(ProjectCont1DSource  ProjectCont1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont1D ProjectCont1DSource)
    TARGET_LINK_LIBRARIES(ProjectCont1D
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g    
    optimized tinyxml debug tinyxml
    )

SET(HelmholtzCont1DSource  Helmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz1D HelmholtzCont1DSource)
    TARGET_LINK_LIBRARIES(Helmholtz1D
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g    
    optimized tinyxml debug tinyxml
    )


#SET(ProjectLoc2DSource  ProjectLoc2D.cpp)
#ADD_NEKTAR_EXECUTABLE(ProjectLoc2D ProjectLoc2DSource)
#    TARGET_LINK_LIBRARIES(ProjectLoc2D
#    optimized MultiRegions debug MultiRegions-g    
#    optimized LocalRegions debug LocalRegions-g    
#    optimized SpatialDomains debug SpatialDomains-g    
#    optimized StdRegions debug StdRegions-g    
#    optimized tinyxml debug tinyxml
#    )

#SET(ProjectCont2DSource  ProjectCont2D.cpp)
#ADD_NEKTAR_EXECUTABLE(ProjectCont2D ProjectCont2DSource)
#    TARGET_LINK_LIBRARIES(ProjectCont2D
#    optimized MultiRegions debug MultiRegions-g    
#    optimized LocalRegions debug LocalRegions-g    
#    optimized SpatialDomains debug SpatialDomains-g    
#    optimized StdRegions debug StdRegions-g    
#    optimized tinyxml debug tinyxml
#    )

IF( USE_FRAMEWORK )
    TARGET_LINK_LIBRARIES(ProjectLoc1D ${FRAMEWORK_LINK_FLAGS})
    TARGET_LINK_LIBRARIES(ProjectCont1D ${FRAMEWORK_LINK_FLAGS})
    TARGET_LINK_LIBRARIES(Helmholtz1D ${FRAMEWORK_LINK_FLAGS})
#    TARGET_LINK_LIBRARIES(ProjectLoc2D ${FRAMEWORK_LINK_FLAGS})
#    TARGET_LINK_LIBRARIES(ProjectCont2D ${FRAMEWORK_LINK_FLAGS})
ENDIF ( USE_FRAMEWORK )

