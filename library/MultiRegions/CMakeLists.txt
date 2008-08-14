SET(MULTI_REGIONS_SOURCES
ContExpList1D.cpp
ContExpList2D.cpp
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
LocalToGlobalMap.cpp
LocalToGlobalDGMap.cpp
LocalToGlobalBaseMap.cpp
LocalToGlobalBndryMap.cpp
LocalToGlobalMap1D.cpp
LocalToGlobalMap2D.cpp
#LocalToGlobalMap3D.cpp
AdvectionContField2D.cpp
)

SET(MULTI_REGIONS_HEADERS
ContExpList1D.h
ContExpList2D.h
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
LocalToGlobalBaseMap.h
LocalToGlobalBndryMap.h
LocalToGlobalDGMap.h
LocalToGlobalMap1D.h
LocalToGlobalMap2D.h
LocalToGlobalMap3D.h
LocalToGlobalMap.h
AdvectionContField2D.h
MultiRegions.hpp
)


ADD_NEKTAR_LIBRARY(MultiRegions STATIC ${MULTI_REGIONS_SOURCES} ${MULTI_REGIONS_HEADERS})

INSTALL(FILES ${MULTI_REGIONS_HEADERS} DESTINATION ${NEKTAR_INCLUDE_DIR}/MultiRegions)





