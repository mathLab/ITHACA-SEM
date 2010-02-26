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

SET(ProjectLoc3DSource  ProjectLoc3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc3D ProjectLoc3DSource)
TARGET_LINK_LIBRARIES(ProjectLoc3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectLoc3D)

SET(ProjectCont3DSource  ProjectCont3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont3D ProjectCont3DSource)
TARGET_LINK_LIBRARIES(ProjectCont3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectCont3D)

SET(ProjectContField3DSource  ProjectContField3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectContField3D ProjectContField3DSource)
TARGET_LINK_LIBRARIES(ProjectContField3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectContField3D)

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
