SET(ProjectLoc1DSource  ProjectLoc1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc1D ProjectLoc1DSource)
    TARGET_LINK_LIBRARIES(ProjectLoc1D
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g    
    optimized tinyxml debug tinyxml
    )
    SET_LAPACK_LINK_LIBRARIES(ProjectLoc1D)

SET(ProjectCont1DSource  ProjectCont1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont1D ProjectCont1DSource)
    TARGET_LINK_LIBRARIES(ProjectCont1D
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g    
    optimized tinyxml debug tinyxml
    )
    SET_LAPACK_LINK_LIBRARIES(ProjectCont1D)

SET(HelmholtzCont1DSource  Helmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz1D HelmholtzCont1DSource)
    TARGET_LINK_LIBRARIES(Helmholtz1D
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g    
    optimized tinyxml debug tinyxml
    )
    SET_LAPACK_LINK_LIBRARIES(Helmholtz1D)

SET(HelmholtzCont2DSource  Helmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz2D HelmholtzCont2DSource)
    TARGET_LINK_LIBRARIES(Helmholtz2D
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g    
    optimized tinyxml debug tinyxml
   )
    SET_LAPACK_LINK_LIBRARIES(Helmholtz2D)

SET(UDGHelmholtz1DSource  UDGHelmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(UDGHelmholtz1D UDGHelmholtz1DSource)
    TARGET_LINK_LIBRARIES(UDGHelmholtz1D
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g    
    optimized tinyxml debug tinyxml
    )
    SET_LAPACK_LINK_LIBRARIES(UDGHelmholtz1D)

#SET(ProjectLoc2DSource  ProjectLoc2D.cpp)
#ADD_NEKTAR_EXECUTABLE(ProjectLoc2D ProjectLoc2DSource)
#    TARGET_LINK_LIBRARIES(ProjectLoc2D
#    optimized MultiRegions debug MultiRegions-g    
#    optimized LocalRegions debug LocalRegions-g    
#    optimized SpatialDomains debug SpatialDomains-g    
#    optimized StdRegions debug StdRegions-g    
#    optimized tinyxml debug tinyxml
#    )
#    SET_LAPACK_LINK_LIBRARIES(ProjectLoc2D)

#SET(ProjectCont2DSource  ProjectCont2D.cpp)
#ADD_NEKTAR_EXECUTABLE(ProjectCont2D ProjectCont2DSource)
#    TARGET_LINK_LIBRARIES(ProjectCont2D
#    optimized MultiRegions debug MultiRegions-g    
#    optimized LocalRegions debug LocalRegions-g    
#    optimized SpatialDomains debug SpatialDomains-g    
#    optimized StdRegions debug StdRegions-g    
#    optimized tinyxml debug tinyxml
#    )
#    SET_LAPACK_LINK_LIBRARIES(ProjectCont2D) 



