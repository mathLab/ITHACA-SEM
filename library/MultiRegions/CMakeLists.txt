SET(MULTI_REGIONS_SOURCES
ContField1D.cpp
ContField2D.cpp
ContField3D.cpp
ContField3DHomogeneous1D.cpp
DisContField1D.cpp
DisContField2D.cpp
DisContField3D.cpp
DisContField3DHomogeneous1D.cpp
ExpList.cpp
ExpList1D.cpp
ExpList2D.cpp
ExpList3D.cpp
ExpListHomogeneous1D.cpp
ExpList2DHomogeneous1D.cpp
ExpList3DHomogeneous1D.cpp
GlobalLinSys.cpp
GlobalLinSysKey.cpp
GlobalLinSysDirect.cpp
GlobalLinSysDirectFull.cpp
GlobalLinSysDirectStaticCond.cpp
GlobalLinSysIterative.cpp
GlobalLinSysIterativeCG.cpp
GlobalMatrix.cpp
GlobalMatrixKey.cpp
LocalMatrixSystem.cpp
LocalToGlobalBaseMap.cpp
LocalToGlobalDGMap.cpp
LocalToGlobalC0ContMap.cpp
GlobalOptimizationParameters.cpp
SubStructuredGraph.cpp
)

SET(MULTI_REGIONS_HEADERS
ContField1D.h
ContField2D.h
ContField3D.h
ContField3DHomogeneous1D.h
DisContField1D.h
DisContField2D.h
DisContField3D.h
DisContField3DHomogeneous1D.h
ExpList.h
ExpList1D.h
ExpList2D.h
ExpList3D.h
ExpListHomogeneous1D.h
ExpList2DHomogeneous1D.h
ExpList3DHomogeneous1D.h
GlobalLinSys.h
GlobalLinSysKey.h
GlobalLinSysDirect.h
GlobalLinSysDirectFull.h
GlobalLinSysIterative.h
GlobalLinSysIterativeCG.h
GlobalMatrix.h
GlobalMatrixKey.h
LocalMatrixSystem.h
LocalToGlobalBaseMap.h
LocalToGlobalDGMap.h
LocalToGlobalC0ContMap.h
MultiRegions.hpp
GlobalOptimizationParameters.h
SubStructuredGraph.h
)


ADD_NEKTAR_LIBRARY(MultiRegions ${NEKTAR_LIBRARY_TYPE} ${MULTI_REGIONS_SOURCES} ${MULTI_REGIONS_HEADERS})
TARGET_LINK_LIBRARIES(MultiRegions
        optimized LibUtilities   debug LibUtilities-g
        optimized SpatialDomains debug SpatialDomains-g
        optimized LocalRegions   debug LocalRegions-g
        optimized StdRegions     debug StdRegions-g
        optimized tinyxml debug tinyxml
	)


INSTALL(FILES ${MULTI_REGIONS_HEADERS} DESTINATION ${NEKTAR_INCLUDE_DIR}/MultiRegions)





