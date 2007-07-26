
SET(GraphSources 
    GraphExample.cpp)
    
SET(MemoryManagerSources 
    MemoryManager.cpp)
    
SET(VecMatSources
    VecMatExample.cpp)

#SET(FoundationSources
#    FoundationDemo.cpp)	
	
SET(NodalTriFeketeSources
	NodalTriFeketeDemo.cpp)	
	
	   
ADD_NEKTAR_EXECUTABLE(Graph GraphSources )
#ADD_NEKTAR_EXECUTABLE(MemoryManager MemoryManagerSources )
ADD_NEKTAR_EXECUTABLE(Matrix VecMatSources )
#ADD_NEKTAR_EXECUTABLE(Foundation FoundationSources )
ADD_NEKTAR_EXECUTABLE(NodalTriFeketeDemo NodalTriFeketeSources )
SET_LAPACK_LINK_LIBRARIES(NodalTriFeketeDemo)

