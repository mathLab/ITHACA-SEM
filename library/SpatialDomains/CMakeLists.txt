SET(SPATIAL_DOMAINS_SOURCES
BoundaryConditions.cpp
#EdgeComponent.cpp
Equation.cpp
GeomFactors.cpp
Geometry1D.cpp
Geometry2D.cpp
Geometry3D.cpp
Geometry.cpp
HexGeom.cpp
InterfaceComponent.cpp
MeshComponents.cpp
MeshGraph.cpp
MeshGraph1D.cpp
MeshGraph2D.cpp
MeshGraph3D.cpp
pchSpatialDomains.cpp
PrismGeom.cpp
PyrGeom.cpp
#QuadFaceComponent.cpp
QuadGeom.cpp
SegGeom.cpp
TetGeom.cpp
#TriFaceComponent.cpp
TriGeom.cpp
)

SET(SPATIAL_DOMAINS_HEADERS
BoundaryConditions.h
Curve.hpp
Equation.h
#EdgeComponent.h
GeomFactors.h
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
ParseUtils.hpp
pchSpatialDomains.h
PrismGeom.h
PyrGeom.h
#QuadFaceComponent.h
QuadGeom.h
SegGeom.h
TetGeom.h
#TriFaceComponent.h
TriGeom.h
SpatialDomains.hpp
)

SET(TINYXML_SOURCES
${TINYXML_SRC_DIR}/tinyxml.cpp
${TINYXML_SRC_DIR}/tinyxmlparser.cpp 
#${TINYXML_SRC_DIR}/xmltest.cpp
${TINYXML_SRC_DIR}/tinyxmlerror.cpp 
${TINYXML_SRC_DIR}/tinystr.cpp
)

SET(TINYXML_HEADERS
${TINYXML_INCLUDE_DIR}/tinystr.h
${TINYXML_INCLUDE_DIR}/tinyxml.h
)


ADD_NEKTAR_LIBRARY(SpatialDomains STATIC ${SPATIAL_DOMAINS_SOURCES} ${SPATIAL_DOMAINS_HEADERS})

INSTALL(FILES ${SPATIAL_DOMAINS_HEADERS} DESTINATION ${NEKTAR_INCLUDE_DIR}/SpatialDomains)


ADD_NEKTAR_LIBRARY(tinyxml STATIC ${TINYXML_SOURCES} ${TINYXML_HEADERS})
INSTALL(TARGETS tinyxml RUNTIME DESTINATION ${NEKTAR_BIN_DIR} OPTIONAL
		        ARCHIVE DESTINATION ${NEKTAR_LIB_DIR} OPTIONAL
			LIBRARY DESTINATION ${NEKTAR_LIB_DIR} OPTIONAL) 

