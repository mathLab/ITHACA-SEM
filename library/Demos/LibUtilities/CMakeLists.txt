
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

