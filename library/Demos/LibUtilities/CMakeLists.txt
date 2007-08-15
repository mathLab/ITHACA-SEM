
SET(GraphSources 
    GraphExample.cpp)
    
SET(MemoryManagerSources 
    MemoryManager.cpp)
    
SET(VecMatSources
    VecMatExample.cpp)

SET(FoundationSources
    FoundationDemo.cpp)	
	
SET(NodalTriFeketeSources
	NodalTriFeketeDemo.cpp)	

SET(NodalTriElecSources
    NodalTriElecDemo.cpp)	

SET(NodalTetElecSources
	NodalTetElecDemo.cpp)			

	   
ADD_NEKTAR_EXECUTABLE(Graph GraphSources )
#ADD_NEKTAR_EXECUTABLE(MemoryManager MemoryManagerSources )
ADD_NEKTAR_EXECUTABLE(Matrix VecMatSources )
ADD_NEKTAR_EXECUTABLE(FoundationDemo FoundationSources )
ADD_NEKTAR_EXECUTABLE(NodalTriFeketeDemo NodalTriFeketeSources )
ADD_NEKTAR_EXECUTABLE(NodalTriElecDemo NodalTriElecSources)
ADD_NEKTAR_EXECUTABLE(NodalTetElecDemo NodalTetElecSources)


