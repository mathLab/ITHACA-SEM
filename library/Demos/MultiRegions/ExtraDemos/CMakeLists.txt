ADD_NEKTAR_EXECUTABLE(ProjectLoc1D
    COMPONENT demos-extra DEPENDS MultiRegions SOURCES ProjectLoc1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc2D
    COMPONENT demos-extra DEPENDS MultiRegions SOURCES ProjectLoc2D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectLoc3D
    COMPONENT demos-extra DEPENDS MultiRegions SOURCES ProjectLoc3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont1D
    COMPONENT demos-extra DEPENDS MultiRegions SOURCES ProjectCont1D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont2D
    COMPONENT demos-extra DEPENDS MultiRegions SOURCES ProjectCont2D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectCont3D
    COMPONENT demos-extra DEPENDS MultiRegions SOURCES ProjectCont3D.cpp)
ADD_NEKTAR_EXECUTABLE(ProjectContField3D
    COMPONENT demos-extra DEPENDS MultiRegions SOURCES ProjectContField3D.cpp)
ADD_NEKTAR_EXECUTABLE(SteadyLinearAdvectionReaction2D
    COMPONENT demos-extra DEPENDS MultiRegions SOURCES SteadyLinearAdvectionReaction2D.cpp)
ADD_NEKTAR_EXECUTABLE(PostProcessing
    COMPONENT demos-extra DEPENDS MultiRegions SOURCES PostProcessing.cpp)

ADD_NEKTAR_TEST(ProjectionCont2D_P7)
IF (NEKTAR_USE_MPI)
    ADD_NEKTAR_TEST(ProjectionCont2D_P7_par)
ENDIF (NEKTAR_USE_MPI)
