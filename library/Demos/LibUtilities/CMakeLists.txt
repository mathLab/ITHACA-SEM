
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

SET(NodalTriEvenlySpacedSources
    NodalTriEvenlySpacedDemo.cpp)

SET(NodalTetEvenlySpacedSources
    NodalTetEvenlySpacedDemo.cpp)

SET(NodalTetElecSources
	NodalTetElecDemo.cpp)
   

   		
ADD_NEKTAR_EXECUTABLE(Graph GraphSources )
#ADD_NEKTAR_EXECUTABLE(MemoryManager MemoryManagerSources )
ADD_NEKTAR_EXECUTABLE(Matrix VecMatSources)

ADD_NEKTAR_EXECUTABLE(FoundationDemo FoundationSources )
SET_LAPACK_LINK_LIBRARIES(FoundationDemo)

ADD_NEKTAR_EXECUTABLE(NodalTriFeketeDemo NodalTriFeketeSources )
SET_LAPACK_LINK_LIBRARIES(NodalTriFeketeDemo)

ADD_NEKTAR_EXECUTABLE(NodalTriElecDemo NodalTriElecSources)
SET_LAPACK_LINK_LIBRARIES(NodalTriElecDemo)

ADD_NEKTAR_EXECUTABLE(NodalTriEvenlySpacedDemo NodalTriEvenlySpacedSources)
SET_LAPACK_LINK_LIBRARIES(NodalTriEvenlySpacedDemo)

ADD_NEKTAR_EXECUTABLE(NodalTetEvenlySpacedDemo NodalTetEvenlySpacedSources)
SET_LAPACK_LINK_LIBRARIES(NodalTetEvenlySpacedDemo)

ADD_NEKTAR_EXECUTABLE(NodalTetElecDemo NodalTetElecSources)
SET_LAPACK_LINK_LIBRARIES(NodalTetElecDemo)



