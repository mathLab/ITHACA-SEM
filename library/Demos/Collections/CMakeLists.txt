SET(LinkLibraries MultiRegions)

SET(CollectionTest2DSource  CollectionTest2D.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTest2D demos CollectionTest2DSource)
TARGET_LINK_LIBRARIES(CollectionTest2D ${LinkLibraries})

SET(CollectionTest3DSource  CollectionTest3D.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTest3D demos CollectionTest3DSource)
TARGET_LINK_LIBRARIES(CollectionTest3D ${LinkLibraries})

#ADD_NEKTAR_TEST(LinearAdvDiffReact2D_P7_Modes)

IF (NEKTAR_USE_MPI)
#    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_xxt_full)
ENDIF (NEKTAR_USE_MPI)
