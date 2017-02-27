
#SET(GraphSources
#    GraphExample.cpp)

SET(MemoryManagerSources
    MemoryManager.cpp)

SET(PartitionAnalyseSources
    PartitionAnalyse.cpp)

SET(FoundationSources
    FoundationDemo.cpp)

SET(NodalDemoSources NodalDemo.cpp)

SET(TimeIntegrationDemoSources
	TimeIntegrationDemo.cpp)

SET(FieldIOBenchmarkerSources
	FieldIOBenchmarker.cpp)

#ADD_NEKTAR_EXECUTABLE(Graph demos GraphSources )
#SET_LAPACK_LINK_LIBRARIES(Graph)

#ADD_NEKTAR_EXECUTABLE(MemoryManager MemoryManagerSources )

ADD_NEKTAR_EXECUTABLE(PartitionAnalyse demos PartitionAnalyseSources)
TARGET_LINK_LIBRARIES(PartitionAnalyse LibUtilities)

ADD_NEKTAR_EXECUTABLE(FoundationDemo demos FoundationSources )
TARGET_LINK_LIBRARIES(FoundationDemo LibUtilities)

ADD_NEKTAR_EXECUTABLE(NodalDemo demos NodalDemoSources)
TARGET_LINK_LIBRARIES(NodalDemo LibUtilities)

ADD_NEKTAR_EXECUTABLE(TimeIntegrationDemo demos TimeIntegrationDemoSources)
TARGET_LINK_LIBRARIES(TimeIntegrationDemo LibUtilities)

ADD_NEKTAR_TEST(NodalDemo_Tri_Deriv_P8)
ADD_NEKTAR_TEST(NodalDemo_Tri_Integral_P6)
ADD_NEKTAR_TEST(NodalDemo_Tri_Interp_P7)
ADD_NEKTAR_TEST(NodalDemo_Prism_Deriv_P8)
ADD_NEKTAR_TEST(NodalDemo_Prism_Integral_P6)
ADD_NEKTAR_TEST(NodalDemo_Prism_Interp_P7)
ADD_NEKTAR_TEST(NodalDemo_Tet_Deriv_P8)
ADD_NEKTAR_TEST(NodalDemo_Tet_Integral_P6)
ADD_NEKTAR_TEST(NodalDemo_Tet_Interp_P7)

IF(NEKTAR_USE_HDF5)
    ADD_NEKTAR_EXECUTABLE(FieldIOBenchmarker demos FieldIOBenchmarkerSources)
    TARGET_LINK_LIBRARIES(FieldIOBenchmarker LibUtilities)
ENDIF()
