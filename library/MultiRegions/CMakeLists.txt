SET(MULTI_REGIONS_SOURCES
ContExpList1D.cpp
ContExpList3D.cpp
ContField1D.cpp
ContField2D.cpp
ContField3D.cpp
DisContField1D.cpp
DisContField2D.cpp
ExpList1D.cpp
ExpList2D.cpp
ExpList3D.cpp
ExpList.cpp
GenExpList1D.cpp
GlobalLinSys.cpp
GlobalLinSysKey.cpp
GlobalMatrix.cpp
GlobalMatrixKey.cpp
LocalToGlobalBaseMap.cpp
LocalToGlobalDGMap.cpp
LocalToGlobalC0ContMap.cpp
GlobalOptimizationParameters.cpp
)

SET(MULTI_REGIONS_HEADERS
ContExpList1D.h
ContExpList3D.h
ContField1D.h
ContField2D.h
ContField3D.h
DisContField1D.h
DisContField2D.h
ExpList1D.h
ExpList2D.h
ExpList3D.h
ExpList.h
GenExpList1D.h
GlobalLinSys.h
GlobalLinSysKey.h
GlobalMatrix.h
GlobalMatrixKey.h
LocalToGlobalBaseMap.h
LocalToGlobalDGMap.h
LocalToGlobalC0ContMap.h
#LocalToGlobalBndryMap.h
#LocalToGlobalMap1D.h
#LocalToGlobalMap2D.h
#LocalToGlobalMap3D.h
#LocalToGlobalMap.h
MultiRegions.hpp
GlobalOptimizationParameters.h
)


ADD_NEKTAR_LIBRARY(MultiRegions STATIC ${MULTI_REGIONS_SOURCES} ${MULTI_REGIONS_HEADERS})

INSTALL(FILES ${MULTI_REGIONS_HEADERS} DESTINATION ${NEKTAR_INCLUDE_DIR}/MultiRegions)





