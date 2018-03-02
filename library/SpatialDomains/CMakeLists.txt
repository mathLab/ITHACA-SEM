SET(SPATIAL_DOMAINS_SOURCES
Conditions.cpp
GeomFactors.cpp
Geometry0D.cpp
Geometry1D.cpp
Geometry2D.cpp
Geometry3D.cpp
Geometry.cpp
HexGeom.cpp
MeshGraph.cpp
MeshGraphXml.cpp
MeshGraphXmlCompressed.cpp
MeshPartition.cpp
MeshPartitionMetis.cpp
PrismGeom.cpp
PyrGeom.cpp
QuadGeom.cpp
SegGeom.cpp
PointGeom.cpp
TetGeom.cpp
TriGeom.cpp
)

SET(SPATIAL_DOMAINS_HEADERS
Conditions.h
Curve.hpp
GeomFactors.h
Geometry0D.h
Geometry1D.h
Geometry2D.h
Geometry3D.h
Geometry.h
HexGeom.h
MeshEntities.hpp
MeshGraph.h
MeshGraphXml.h
MeshGraphXmlCompressed.h
MeshPartition.h
MeshPartitionMetis.h
PrismGeom.h
PyrGeom.h
QuadGeom.h
SegGeom.h
PointGeom.h
SpatialDomainsDeclspec.h
TetGeom.h
TriGeom.h
SpatialDomains.hpp
)

IF(NEKTAR_USE_HDF5)
	SET(SPATIAL_DOMAINS_HEADERS ${SPATIAL_DOMAINS_HEADERS} MeshGraphHDF5.h)
	SET(SPATIAL_DOMAINS_SOURCES ${SPATIAL_DOMAINS_SOURCES} MeshGraphHDF5.cpp)
ENDIF()

IF(NEKTAR_USE_SCOTCH)
	SET(SPATIAL_DOMAINS_HEADERS ${SPATIAL_DOMAINS_HEADERS} MeshPartitionScotch.h)
	SET(SPATIAL_DOMAINS_SOURCES ${SPATIAL_DOMAINS_SOURCES} MeshPartitionScotch.cpp)
ENDIF()

ADD_DEFINITIONS(-DSPATIAL_DOMAINS_EXPORTS)

ADD_NEKTAR_LIBRARY(SpatialDomains
    SOURCES ${SPATIAL_DOMAINS_SOURCES}
    HEADERS ${SPATIAL_DOMAINS_HEADERS}
    DEPENDS StdRegions
    SUMMARY "Nektar++ SpatialDomains library"
    DESCRIPTION "This library provides the mappings between reference regions and physical regions in the domain.")

ADD_DEPENDENCIES(SpatialDomains modmetis-5.1.0)

# Metis
TARGET_LINK_LIBRARIES(SpatialDomains LINK_PRIVATE ${METIS_LIB})

IF( NEKTAR_USE_SCOTCH )
    TARGET_LINK_LIBRARIES(SpatialDomains LINK_PRIVATE
        ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY})
    ADD_DEPENDENCIES(SpatialDomains scotch-6.0.0)
ENDIF ()

IF (NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/SpatialDomains
    COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
