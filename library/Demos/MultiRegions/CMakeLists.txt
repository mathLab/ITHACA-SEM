SET(LinkLibraries 
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g 
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
    optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
    optimized ${Boost_ZLIB_LIBRARY_RELEASE} debug ${Boost_ZLIB_LIBRARY_DEBUG}
    optimized tinyxml debug tinyxml
    optimized sbtk    debug sbtk
)

SET(FldToTec2DSource  FldToTec2D.cpp)
ADD_NEKTAR_EXECUTABLE(FldToTec2D FldToTec2DSource)
TARGET_LINK_LIBRARIES(FldToTec2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(FldToTec2D)

SET(FldToGmsh2DSource  FldToGmsh.cpp)
ADD_NEKTAR_EXECUTABLE(FldToGmsh2D FldToGmsh2DSource)
TARGET_LINK_LIBRARIES(FldToGmsh2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(FldToGmsh2D)

SET(DatToFldSource  DatToFld.cpp)
ADD_NEKTAR_EXECUTABLE(DatToFld DatToFldSource)
TARGET_LINK_LIBRARIES(DatToFld ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(DatToFld)

SET(ProjectLoc1DSource  ProjectLoc1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc1D ProjectLoc1DSource)
TARGET_LINK_LIBRARIES(ProjectLoc1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectLoc1D)


SET(ProjectCont1DSource  ProjectCont1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont1D ProjectCont1DSource)
TARGET_LINK_LIBRARIES(ProjectCont1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectCont1D)


SET(HelmholtzCont1DSource  Helmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz1D HelmholtzCont1DSource)
TARGET_LINK_LIBRARIES(Helmholtz1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz1D)


SET(HelmholtzCont2DSource  Helmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz2D HelmholtzCont2DSource)
TARGET_LINK_LIBRARIES(Helmholtz2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz2D)

SET(HelmholtzCont3DSource  Helmholtz3D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3D HelmholtzCont3DSource)
TARGET_LINK_LIBRARIES(Helmholtz3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz3D)


SET(HDGHelmholtz1DSource  HDGHelmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz1D HDGHelmholtz1DSource)
TARGET_LINK_LIBRARIES(HDGHelmholtz1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz1D)


SET(HDGHelmholtz2DSource  HDGHelmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz2D HDGHelmholtz2DSource)
TARGET_LINK_LIBRARIES(HDGHelmholtz2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz2D)

# SET(HDGHelmholtz3DSource  HDGHelmholtz3D.cpp)
# ADD_NEKTAR_EXECUTABLE(HDGHelmholtz3D HDGHelmholtz3DSource)
# TARGET_LINK_LIBRARIES(HDGHelmholtz3D ${LinkLibraries})
# SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz3D)

SET(ProjectLoc2DSource  ProjectLoc2D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc2D ProjectLoc2DSource)
TARGET_LINK_LIBRARIES(ProjectLoc2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectLoc2D)

SET(ProjectCont2DSource  ProjectCont2D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont2D ProjectCont2DSource)
TARGET_LINK_LIBRARIES(ProjectCont2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectCont2D) 


SET(AdvectionSource  Advection.cpp)
ADD_NEKTAR_EXECUTABLE(Advection AdvectionSource)
TARGET_LINK_LIBRARIES(Advection ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Advection)

SET(LinearAdvectionCont2DSource  LinearAdvection2D.cpp)
ADD_NEKTAR_EXECUTABLE(LinearAdvection2D LinearAdvectionCont2DSource)
TARGET_LINK_LIBRARIES(LinearAdvection2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LinearAdvection2D)

SET(EigValsLinearAdvectionCont2DSource  EigValsLinearAdvection2D.cpp)
ADD_NEKTAR_EXECUTABLE(EigValsLinearAdvection2D EigValsLinearAdvectionCont2DSource)
TARGET_LINK_LIBRARIES(EigValsLinearAdvection2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(EigValsLinearAdvection2D)

SET(PostProcessingSource  PostProcessing.cpp)
ADD_NEKTAR_EXECUTABLE(PostProcessing PostProcessingSource)
TARGET_LINK_LIBRARIES(PostProcessing ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(PostProcessing)
