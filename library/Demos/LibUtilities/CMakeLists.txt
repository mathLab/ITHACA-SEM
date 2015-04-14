
#SET(GraphSources 
#    GraphExample.cpp)
    
SET(MemoryManagerSources 
    MemoryManager.cpp)
    
SET(PartitionAnalyseSources
    PartitionAnalyse.cpp)
    
SET(FoundationSources
    FoundationDemo.cpp)	
	
SET(NodalTriFeketeSources
	NodalTriFeketeDemo.cpp)

SET(NodalTriElecSources
    NodalTriElecDemo.cpp)

SET(NodalTriEvenlySpacedSources
    NodalTriEvenlySpacedDemo.cpp)

SET(NodalTetEvenlySpacedSources
    NodalTetEvenlySpacedDemo.cpp)

SET(NodalTetElecSources
	NodalTetElecDemo.cpp)

SET(TimeIntegrationDemoSources
	TimeIntegrationDemo.cpp)


#ADD_NEKTAR_EXECUTABLE(Graph demos GraphSources )
#SET_LAPACK_LINK_LIBRARIES(Graph)

#ADD_NEKTAR_EXECUTABLE(MemoryManager MemoryManagerSources )

ADD_NEKTAR_EXECUTABLE(PartitionAnalyse demos PartitionAnalyseSources)
TARGET_LINK_LIBRARIES(PartitionAnalyse LibUtilities)

ADD_NEKTAR_EXECUTABLE(FoundationDemo demos FoundationSources )
TARGET_LINK_LIBRARIES(FoundationDemo LibUtilities)

ADD_NEKTAR_EXECUTABLE(NodalTriFeketeDemo demos NodalTriFeketeSources )
TARGET_LINK_LIBRARIES(NodalTriFeketeDemo LibUtilities)

ADD_NEKTAR_EXECUTABLE(NodalTriElecDemo demos NodalTriElecSources)
TARGET_LINK_LIBRARIES(NodalTriElecDemo LibUtilities)

ADD_NEKTAR_EXECUTABLE(NodalTriEvenlySpacedDemo demos NodalTriEvenlySpacedSources)
TARGET_LINK_LIBRARIES(NodalTriEvenlySpacedDemo LibUtilities)

ADD_NEKTAR_EXECUTABLE(NodalTetEvenlySpacedDemo demos NodalTetEvenlySpacedSources)
TARGET_LINK_LIBRARIES(NodalTetEvenlySpacedDemo LibUtilities)

ADD_NEKTAR_EXECUTABLE(NodalTetElecDemo demos NodalTetElecSources)
TARGET_LINK_LIBRARIES(NodalTetElecDemo LibUtilities)

ADD_NEKTAR_EXECUTABLE(TimeIntegrationDemo demos TimeIntegrationDemoSources)
TARGET_LINK_LIBRARIES(TimeIntegrationDemo LibUtilities)

