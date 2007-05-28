SET(MULTI_REGIONS_SOURCES
ContExpList1D.cpp
#DisContExpList1D.cpp
ExpList1D.cpp
#ExpList2D.cpp
#ExpList3D.cpp
ExpList.cpp
LocalToGlobalMap.cpp
LocalToGlobalMap1D.cpp
#LocalToGlobalMap2D.cpp
#LocalToGlobalMap3D.cpp
)

SET(MULTI_REGIONS_HEADERS
ContExpList1D.h
DisContExpList1D.h
ExpList1D.h
ExpList2D.h
ExpList3D.h
ExpList.h
LocalToGlobalMap1D.h
LocalToGlobalMap2D.h
LocalToGlobalMap3D.h
LocalToGlobalMap.h
MultiRegions.hpp
)


ADD_NEKTAR_LIBRARY(MultiRegions STATIC ${MULTI_REGIONS_SOURCES} ${MULTI_REGIONS_HEADERS})

INSTALL_FILES(${NEKTAR_INCLUDE_DIR}/MultiRegions FILES ${MULTI_REGIONS_HEADERS})




