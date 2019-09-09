SET(LOCAL_REGIONS_SOURCES
Expansion.cpp
Expansion0D.cpp
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
Expansion0D.h
Expansion1D.h
Expansion2D.h
Expansion3D.h
HexExp.h
LinSys.hpp
#LocalRegions.h
LocalRegions.hpp
LocalRegionsDeclspec.h
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

ADD_DEFINITIONS(-DLOCAL_REGIONS_EXPORTS)

ADD_NEKTAR_LIBRARY(LocalRegions
    SOURCES ${LOCAL_REGIONS_SOURCES}
    HEADERS ${LOCAL_REGIONS_HEADERS}
    DEPENDS SpatialDomains
    SUMMARY "Nektar++ LocalRegions library"
    DESCRIPTION "This library provides physical space expansions on the various supported regions.")

IF (NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/LocalRegions
    COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
