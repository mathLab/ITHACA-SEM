SET(LOCAL_REGIONS_SOURCES
Expansion.cpp
Expansion1D.cpp
Expansion2D.cpp
Expansion3D.cpp
QuadExp.cpp
HexExp.cpp
MatrixKey.cpp
NodalTetExp.cpp
NodalTriExp.cpp
PointExp.cpp
PrismExp.cpp
PyrExp.cpp
SegExp.cpp
TetExp.cpp
TriExp.cpp
)

SET(LOCAL_REGIONS_HEADERS
Expansion.h
Expansion1D.h
Expansion2D.h
Expansion3D.h
HexExp.h
LinSys.hpp
LocalRegions.h
LocalRegions.hpp
MatrixKey.h
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
ADD_NEKTAR_LIBRARY(LocalRegions ${NEKTAR_LIBRARY_TYPE} ${LOCAL_REGIONS_SOURCES} ${LOCAL_REGIONS_HEADERS})
SETUP_PRECOMPILED_HEADERS(LOCAL_REGIONS_SOURCES PrecompiledHeaderName)

INSTALL(FILES ${LOCAL_REGIONS_HEADERS} DESTINATION ${NEKTAR_INCLUDE_DIR}/LocalRegions)
