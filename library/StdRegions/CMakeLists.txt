SET(STD_REGIONS_SOURCES
BasisManager.cpp
NodalBasisManager.cpp
PolyManager.cpp
StdBasis.cpp
StdExpansion1D.cpp
StdExpansion2D.cpp
StdExpansion3D.cpp
StdExpansion.cpp
StdExpMap.cpp
StdHexExp.cpp
StdMatrix.cpp
StdNodalTetExp.cpp
StdNodalTriExp.cpp
StdPrismExp.cpp
StdPyrExp.cpp
StdQuadExp.cpp
StdSegExp.cpp
StdTetExp.cpp
StdTriExp.cpp
)

SET(STD_REGIONS_HEADERS
BasisManager.h
NodalBasisManager.h
PolyManager.h
SpatialDomainsDeclarations.h
StdBasis.h
StdExpansion1D.h
StdExpansion2D.h
StdExpansion3D.h
StdExpansion.h
StdExpMap.h
StdHexExp.h
StdMatrix.h
StdNodalTetExp.h
StdNodalTriExp.h
StdPoints.h
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

