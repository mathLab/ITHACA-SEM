SET(LinkLibraries MultiRegions)

SET(CollectionTimingSource  CollectionTiming.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTiming demos CollectionTimingSource)
TARGET_LINK_LIBRARIES(CollectionTiming ${LinkLibraries})

#ADD_NEKTAR_TEST(LinearAdvDiffReact2D_P7_Modes)

IF (NEKTAR_USE_MPI)
#    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_xxt_full)
ENDIF (NEKTAR_USE_MPI)
