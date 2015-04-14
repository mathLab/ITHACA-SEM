SET(LinkLibraries MultiRegions)

SET(CollectionTest2DSource  CollectionTest2D.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTest2D demos CollectionTest2DSource)
TARGET_LINK_LIBRARIES(CollectionTest2D ${LinkLibraries})

SET(CollectionTiming2DSource  CollectionTiming2D.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTiming2D demos CollectionTiming2DSource)
TARGET_LINK_LIBRARIES(CollectionTiming2D ${LinkLibraries})

SET(CollectionTest3DSource  CollectionTest3D.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTest3D demos CollectionTest3DSource)
TARGET_LINK_LIBRARIES(CollectionTest3D ${LinkLibraries})

SET(CollectionTiming3DSource  CollectionTiming3D.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTiming3D demos CollectionTiming3DSource)
TARGET_LINK_LIBRARIES(CollectionTiming3D ${LinkLibraries})

SET(CollectionTimingSource  CollectionTiming.cpp)
ADD_NEKTAR_EXECUTABLE(CollectionTiming demos CollectionTimingSource)
TARGET_LINK_LIBRARIES(CollectionTiming ${LinkLibraries})

#ADD_NEKTAR_TEST(LinearAdvDiffReact2D_P7_Modes)

IF (NEKTAR_USE_MPI)
#    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_xxt_full)
ENDIF (NEKTAR_USE_MPI)
