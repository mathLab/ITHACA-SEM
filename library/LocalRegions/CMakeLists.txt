SET(LOCAL_REGIONS_SOURCES
GenSegExp.cpp
HexExp.cpp
MatrixKey.cpp
#MetricRelatedInfo.cpp
NodalTetExp.cpp
NodalTriExp.cpp
PointExp.cpp
PrismExp.cpp
PyrExp.cpp
QuadExp.cpp
SegExp.cpp
TetExp.cpp
TriExp.cpp
)

SET(LOCAL_REGIONS_HEADERS
GenSegExp.h
HexExp.h
LinSys.hpp
LocalRegions.h
LocalRegions.hpp
MatrixKey.h
#MetricRelatedInfo.h
NodalTetExp.h
NodalTriExp.h
PointExp.h
PrismExp.h
PyrExp.h
QuadExp.h
SegExp.h
TetExp.h
TriExp.h
)


SET(PrecompiledHeaderName LocalRegions/LocalRegions.h)
ADD_NEKTAR_LIBRARY(LocalRegions STATIC ${LOCAL_REGIONS_SOURCES} ${LOCAL_REGIONS_HEADERS})
SETUP_PRECOMPILED_HEADERS(LOCAL_REGIONS_SOURCES PrecompiledHeaderName)

INSTALL(FILES ${LOCAL_REGIONS_HEADERS} DESTINATION ${NEKTAR_INCLUDE_DIR}/LocalRegions)
