
SET(GraphSources 
    GraphExample.cpp)
    
SET(MemoryManagerSources 
    MemoryManager.cpp)
    
SET(VecMatSources
    VecMatExample.cpp)
    
ADD_NEKTAR_EXECUTABLE(Graph GraphSources )
ADD_NEKTAR_EXECUTABLE(MemoryManager MemoryManagerSources )
ADD_NEKTAR_EXECUTABLE(Matrix VecMatSources )

