SET(SPATIALDOMAINS_SOURCES
    SpatialDomains.cpp
    Geometry.cpp
    Geometry1D.cpp
    Geometry2D.cpp
    MeshGraph.cpp
    SegGeom.cpp
    QuadGeom.cpp
    TriGeom.cpp
    )

ADD_NEKPY_LIBRARY(SpatialDomains SPATIALDOMAINS_SOURCES)
