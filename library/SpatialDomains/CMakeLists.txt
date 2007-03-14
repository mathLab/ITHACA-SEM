SET(SPATIAL_DOMAINS_SOURCES
EdgeComponent.cpp
GeoFac.cpp
Geometry1D.cpp
#Geometry2D.cpp
#Geometry3D.cpp
Geometry.cpp
#HexGeom.cpp
InterfaceComponent.cpp
MeshComponents.cpp
MeshGraph.cpp
MeshGraph1D.cpp
#MeshGraph2D.cpp
#MeshGraph3D.cpp
pchSpatialDomains.cpp
#PrismGeom.cpp
#PyrGeom.cpp
#QuadFaceComponent.cpp
#QuadGeom.cpp
SegGeom.cpp
#TetGeom.cpp
#TriFaceComponent.cpp
#TriGeom.cpp
)

SET(SPATIAL_DOMAINS_HEADERS
EdgeComponent.h
GeoFac.h
Geometry1D.h
Geometry2D.h
Geometry3D.h
Geometry.h
HexGeom.h
InterfaceComponent.h
MeshComponents.h
MeshGraph1D.h
MeshGraph2D.h
MeshGraph3D.h
MeshGraph.h
pchSpatialDomains.h
PrismGeom.h
PyrGeom.h
QuadFaceComponent.h
QuadGeom.h
SegGeom.h
TetGeom.h
TriFaceComponent.h
TriGeom.h
SpatialDomains.hpp
)

SET(TINYXML_SOURCES
${TINYXML_SRC_DIR}/tinyxml.cpp
${TINYXML_SRC_DIR}/tinyxmlparser.cpp 
${TINYXML_SRC_DIR}/xmltest.cpp
${TINYXML_SRC_DIR}/tinyxmlerror.cpp 
${TINYXML_SRC_DIR}/tinystr.cpp
)

SET(TINYXML_HEADERS
${TINYXML_INCLUDE_DIR}/tinystr.h
${TINYXML_INCLUDE_DIR}/tinyxml.h
)


ADD_NEKTAR_LIBRARY(SpatialDomains STATIC ${SPATIAL_DOMAINS_SOURCES} ${SPATIAL_DOMAINS_HEADERS})

INSTALL_FILES(${NEKTAR_INCLUDE_DIR}/SpatialDomains FILES ${SPATIAL_DOMAINS_HEADERS})

ADD_LIBRARY(tinyxml STATIC ${TINYXML_SOURCES} ${TINYXML_HEADERS})


