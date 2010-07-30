SET(MULTI_REGIONS_SOURCES
ContField1D.cpp
ContField2D.cpp
ContField3D.cpp
DisContField1D.cpp
DisContField2D.cpp
DisContField3D.cpp
ExpList1D.cpp
ExpList2D.cpp
ExpList3D.cpp
ExpList3DHomogeneous1D.cpp
ExpList.cpp
GlobalLinSys.cpp
GlobalLinSysKey.cpp
GlobalMatrix.cpp
GlobalMatrixKey.cpp
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
DisContField1D.h
DisContField2D.h
DisContField3D.h
ExpList1D.h
ExpList2D.h
ExpList3D.h
ExpList3DHomogeneous1D.h
ExpList.h
GlobalLinSys.h
GlobalLinSysKey.h
GlobalMatrix.h
GlobalMatrixKey.h
LocalToGlobalBaseMap.h
LocalToGlobalDGMap.h
LocalToGlobalC0ContMap.h
MultiRegions.hpp
GlobalOptimizationParameters.h
SubStructuredGraph.h
)


ADD_NEKTAR_LIBRARY(MultiRegions ${NEKTAR_LIBRARY_TYPE} ${MULTI_REGIONS_SOURCES} ${MULTI_REGIONS_HEADERS})

INSTALL(FILES ${MULTI_REGIONS_HEADERS} DESTINATION ${NEKTAR_INCLUDE_DIR}/MultiRegions)





