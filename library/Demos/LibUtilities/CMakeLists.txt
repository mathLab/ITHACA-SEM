
#SET(GraphSources 
#    GraphExample.cpp)
    
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

SET(TimeIntegrationDemoSources
	TimeIntegrationDemo.cpp)

SET(SvtvpHeaderInlineSources
	Svtvp-header-inline.hpp
	Svtvp-header-inline.cpp)

SET(SvtvpHeaderSources
	Svtvp-header.hpp
	Svtvp-header.cpp)


SET(SvtvpExplicitInstLibSources
	Svtvp-explicit-inst-inc.cpp
	Svtvp-explicit-inst.cpp)
SET(SvtvpExplicitInstLibHeaders
	Svtvp-explicit-inst-inc.hpp
	Svtvp-explicit-inst.hpp)


SET(SvtvpExplicitInstIncSources
    Svtvp-explicit-inst-inc-main.cpp)

SET(SvtvpExplicitInstSources
    Svtvp-explicit-inst-main.cpp)

ADD_NEKTAR_LIBRARY(SvtvpExplicitInstLib lib ${NEKTAR_LIBRARY_TYPE} ${SvtvpExplicitInstLibSources} ${SvtvpExplicitInstLibHeaders})


   		
#ADD_NEKTAR_EXECUTABLE(Graph demos GraphSources )
#SET_LAPACK_LINK_LIBRARIES(Graph)

#ADD_NEKTAR_EXECUTABLE(MemoryManager MemoryManagerSources )
ADD_NEKTAR_EXECUTABLE(Matrix demos VecMatSources)
SET_LAPACK_LINK_LIBRARIES(Matrix)

ADD_NEKTAR_EXECUTABLE(FoundationDemo demos FoundationSources )
SET_LAPACK_LINK_LIBRARIES(FoundationDemo)

ADD_NEKTAR_EXECUTABLE(NodalTriFeketeDemo demos NodalTriFeketeSources )
SET_LAPACK_LINK_LIBRARIES(NodalTriFeketeDemo)

ADD_NEKTAR_EXECUTABLE(NodalTriElecDemo demos NodalTriElecSources)
SET_LAPACK_LINK_LIBRARIES(NodalTriElecDemo)

ADD_NEKTAR_EXECUTABLE(NodalTriEvenlySpacedDemo demos NodalTriEvenlySpacedSources)
SET_LAPACK_LINK_LIBRARIES(NodalTriEvenlySpacedDemo)

ADD_NEKTAR_EXECUTABLE(NodalTetEvenlySpacedDemo demos NodalTetEvenlySpacedSources)
SET_LAPACK_LINK_LIBRARIES(NodalTetEvenlySpacedDemo)

ADD_NEKTAR_EXECUTABLE(NodalTetElecDemo demos NodalTetElecSources)
SET_LAPACK_LINK_LIBRARIES(NodalTetElecDemo)

ADD_NEKTAR_EXECUTABLE(TimeIntegrationDemo demos TimeIntegrationDemoSources)
SET_LAPACK_LINK_LIBRARIES(TimeIntegrationDemo)

ADD_NEKTAR_EXECUTABLE(SvtvpHeaderInlineDemo demos SvtvpHeaderInlineSources)
ADD_NEKTAR_EXECUTABLE(SvtvpHeaderDemo demos SvtvpHeaderSources)

ADD_NEKTAR_EXECUTABLE(SvtvpExplicitInstIncDemo demos SvtvpExplicitInstIncSources)
ADD_NEKTAR_EXECUTABLE(SvtvpExplicitInstDemo    demos SvtvpExplicitInstSources)

TARGET_LINK_LIBRARIES(SvtvpExplicitInstIncDemo SvtvpExplicitInstLib)
TARGET_LINK_LIBRARIES(SvtvpExplicitInstDemo SvtvpExplicitInstLib)

