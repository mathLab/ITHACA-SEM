SET(STD_REGIONS_SOURCES
StdExpansion.cpp
StdExpansion1D.cpp
StdExpansion2D.cpp
StdExpansion3D.cpp
StdHexExp.cpp
#StdNodalTetExp.cpp
StdMatrixKey.cpp
StdNodalTriExp.cpp
StdPrismExp.cpp
StdPyrExp.cpp
StdQuadExp.cpp
StdSegExp.cpp
StdTetExp.cpp
StdTriExp.cpp
OptimizationParametersAccess.cpp
)

SET(STD_REGIONS_HEADERS
SpatialDomainsDeclarations.hpp
StdExpansion1D.h
StdExpansion2D.h
StdExpansion3D.h
StdExpansion.h
StdHexExp.h
StdLinSysKey.hpp
StdMatrixKey.h
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
ElementalOptimizationParameters.hpp
ElementalOptimizationParametersTraits.hpp
OptimizationParametersAccess.h
)


ADD_NEKTAR_LIBRARY(StdRegions ${NEKTAR_LIBRARY_TYPE} ${STD_REGIONS_SOURCES} ${STD_REGIONS_HEADERS})
TARGET_LINK_LIBRARIES(StdRegions
        optimized LibUtilities debug LibUtilities-g
        optimized tinyxml debug tinyxml
	)

INSTALL(FILES ${STD_REGIONS_HEADERS} DESTINATION ${NEKTAR_INCLUDE_DIR}/StdRegions)


#TARGET_LINK_LIBRARIES(StdRegions
#        optimized LibUtilities debug LibUtilities-g
#)

