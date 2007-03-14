SET(LOCAL_REGIONS_SOURCES
#HexExp.cpp
MatrixKey.cpp
#MetricRelatedInfo.cpp
#NodalTetExp.cpp
#NodalTriExp.cpp
#PrismExp.cpp
#PyrExp.cpp
#QuadExp.cpp
SegExp.cpp
#TetExp.cpp
#TriExp.cpp
)

SET(LOCAL_REGIONS_HEADERS
HexExp.h
LinSys.hpp
MatrixKey.h
#MetricRelatedInfo.h
NodalTetExp.h
NodalTriExp.h
PrismExp.h
PyrExp.h
QuadExp.h
SegExp.h
TetExp.h
TriExp.h
)


ADD_NEKTAR_LIBRARY(LocalRegions STATIC ${LOCAL_REGIONS_SOURCES} ${LOCAL_REGIONS_HEADERS})





