SET(LinkLibraries 
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g 
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
    optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
    optimized ${ZLIB_LIBRARY_RELEASE} debug ${ZLIB_LIBRARY_DEBUG}
    optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
    optimized sbtk    debug sbtk
)

SET(FldToTec2DSource  FldToTec2D.cpp)
ADD_NEKTAR_EXECUTABLE(FldToTec2D extra-demos FldToTec2DSource)
TARGET_LINK_LIBRARIES(FldToTec2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(FldToTec2D)

SET(FldToGmsh2DSource  FldToGmsh.cpp)
ADD_NEKTAR_EXECUTABLE(FldToGmsh2D extra-demos FldToGmsh2DSource)
TARGET_LINK_LIBRARIES(FldToGmsh2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(FldToGmsh2D)

SET(DatToFldSource  DatToFld.cpp)
ADD_NEKTAR_EXECUTABLE(DatToFld extra-demos DatToFldSource)
TARGET_LINK_LIBRARIES(DatToFld ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(DatToFld)

SET(ProjectLoc1DSource  ProjectLoc1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc1D extra-demos ProjectLoc1DSource)
TARGET_LINK_LIBRARIES(ProjectLoc1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectLoc1D)

SET(ProjectCont1DSource  ProjectCont1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont1D extra-demos ProjectCont1DSource)
TARGET_LINK_LIBRARIES(ProjectCont1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectCont1D)

SET(ProjectLoc3DSource  ProjectLoc3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc3D extra-demos ProjectLoc3DSource)
TARGET_LINK_LIBRARIES(ProjectLoc3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectLoc3D)

SET(ProjectCont3DSource  ProjectCont3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont3D extra-demos ProjectCont3DSource)
TARGET_LINK_LIBRARIES(ProjectCont3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectCont3D)

SET(ProjectContField3DSource  ProjectContField3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectContField3D extra-demos ProjectContField3DSource)
TARGET_LINK_LIBRARIES(ProjectContField3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectContField3D)

SET(ProjectLoc2DSource  ProjectLoc2D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc2D extra-demos ProjectLoc2DSource)
TARGET_LINK_LIBRARIES(ProjectLoc2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectLoc2D)

SET(ProjectCont2DSource  ProjectCont2D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont2D extra-demos ProjectCont2DSource)
TARGET_LINK_LIBRARIES(ProjectCont2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(ProjectCont2D) 

SET(SteadyLinearAdvectionReactionCont2DSource  SteadyLinearAdvectionReaction2D.cpp)
ADD_NEKTAR_EXECUTABLE(SteadyLinearAdvectionReaction2D extra-demos SteadyLinearAdvectionReactionCont2DSource)
TARGET_LINK_LIBRARIES(SteadyLinearAdvectionReaction2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(SteadyLinearAdvectionReaction2D)

#SET(EigValsLinearAdvectionCont2DSource  EigValsLinearAdvection2D.cpp)
#ADD_NEKTAR_EXECUTABLE(EigValsLinearAdvection2D EigValsLinearAdvectionCont2DSource)
#TARGET_LINK_LIBRARIES(EigValsLinearAdvection2D ${LinkLibraries})
#SET_LAPACK_LINK_LIBRARIES(EigValsLinearAdvection2D)

SET(PostProcessingSource  PostProcessing.cpp)
ADD_NEKTAR_EXECUTABLE(PostProcessing extra-demos PostProcessingSource)
TARGET_LINK_LIBRARIES(PostProcessing ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(PostProcessing)

# Do not build these now as they conflict with library/Timings
#IF( NEKTAR_BUILD_TIMINGS )
#    SET(TimingCGHelmSolve2DSource  TimingCGHelmSolve2D.cpp)
#    ADD_NEKTAR_EXECUTABLE(TimingCGHelmSolve2D TimingCGHelmSolve2DSource)
#    TARGET_LINK_LIBRARIES(TimingCGHelmSolve2D ${LinkLibraries})
#   SET_LAPACK_LINK_LIBRARIES(TimingCGHelmSolve2D)
#    
#    SET(TimingHDGHelmSolve2DSource  TimingHDGHelmSolve2D.cpp)
#    ADD_NEKTAR_EXECUTABLE(TimingHDGHelmSolve2D TimingHDGHelmSolve2DSource)
#    TARGET_LINK_LIBRARIES(TimingHDGHelmSolve2D ${LinkLibraries})
#    SET_LAPACK_LINK_LIBRARIES(TimingHDGHelmSolve2D)
#ENDIF( NEKTAR_BUILD_TIMINGS )
