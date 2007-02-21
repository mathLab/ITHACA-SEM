SET(STD_REGIONS_SOURCES
StdExpansion.cpp
StdExpansion1D.cpp
StdExpManagers.cpp
#StdExpansion2D.cpp
#StdExpansion3D.cpp
StdExpMap.cpp
#StdHexExp.cpp
#StdNodalTetExp.cpp
#StdNodalTriExp.cpp
#StdPrismExp.cpp
#StdPyrExp.cpp
#StdQuadExp.cpp
StdSegExp.cpp
#StdTetExp.cpp
#StdTriExp.cpp
)

SET(STD_REGIONS_HEADERS
SpatialDomainsDeclarations.h
StdExpansion1D.h
StdExpansion2D.h
StdExpansion3D.h
StdExpansion.h
StdExpManagers.h
StdExpMap.h
StdHexExp.h
StdNodalTetExp.h
StdNodalTriExp.h
StdPrismExp.h
StdPyrExp.h
StdQuadExp.h
StdSegExp.h
StdTetExp.h
StdTriExp.h
LocalRegionsDeclarations.hpp
SpatialDomainsDeclarations.hpp
StdRegions.hpp
)


ADD_NEKTAR_LIBRARY(StdRegions STATIC ${STD_REGIONS_SOURCES} ${STD_REGIONS_HEADERS})

INSTALL_FILES(${NEKTAR_INCLUDE_DIR}/StdRegions FILES ${STD_REGIONS_HEADERS})


TARGET_LINK_LIBRARIES(StdRegions
        optimized LibUtilities debug LibUtilities-g
        optimized ${BoostThreadName} debug ${BoostThreadDebugName})

