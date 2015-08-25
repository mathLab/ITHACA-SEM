SET(LinkLibraries MultiRegions)

SET(ProjectLoc1DSource  ProjectLoc1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc1D extra-demos ProjectLoc1DSource)
TARGET_LINK_LIBRARIES(ProjectLoc1D ${LinkLibraries})

SET(ProjectCont1DSource  ProjectCont1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont1D extra-demos ProjectCont1DSource)
TARGET_LINK_LIBRARIES(ProjectCont1D ${LinkLibraries})

SET(ProjectLoc3DSource  ProjectLoc3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc3D extra-demos ProjectLoc3DSource)
TARGET_LINK_LIBRARIES(ProjectLoc3D ${LinkLibraries})

SET(ProjectCont3DSource  ProjectCont3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont3D extra-demos ProjectCont3DSource)
TARGET_LINK_LIBRARIES(ProjectCont3D ${LinkLibraries})

SET(ProjectContField3DSource  ProjectContField3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectContField3D extra-demos ProjectContField3DSource)
TARGET_LINK_LIBRARIES(ProjectContField3D ${LinkLibraries})

SET(ProjectLoc2DSource  ProjectLoc2D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc2D extra-demos ProjectLoc2DSource)
TARGET_LINK_LIBRARIES(ProjectLoc2D ${LinkLibraries})

SET(ProjectCont2DSource  ProjectCont2D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont2D extra-demos ProjectCont2DSource)
TARGET_LINK_LIBRARIES(ProjectCont2D ${LinkLibraries})

SET(SteadyLinearAdvectionReactionCont2DSource  SteadyLinearAdvectionReaction2D.cpp)
ADD_NEKTAR_EXECUTABLE(SteadyLinearAdvectionReaction2D extra-demos SteadyLinearAdvectionReactionCont2DSource)
TARGET_LINK_LIBRARIES(SteadyLinearAdvectionReaction2D ${LinkLibraries})

#SET(EigValsLinearAdvectionCont2DSource  EigValsLinearAdvection2D.cpp)
#ADD_NEKTAR_EXECUTABLE(EigValsLinearAdvection2D EigValsLinearAdvectionCont2DSource)
#TARGET_LINK_LIBRARIES(EigValsLinearAdvection2D ${LinkLibraries})
#SET_LAPACK_LINK_LIBRARIES(EigValsLinearAdvection2D)

SET(PostProcessingSource  PostProcessing.cpp)
ADD_NEKTAR_EXECUTABLE(PostProcessing extra-demos PostProcessingSource)
TARGET_LINK_LIBRARIES(PostProcessing ${LinkLibraries})

ADD_NEKTAR_TEST(ProjectionCont2D_P7)

IF (NEKTAR_USE_MPI)
    ADD_NEKTAR_TEST(ProjectionCont2D_P7_par)
ENDIF (NEKTAR_USE_MPI)
