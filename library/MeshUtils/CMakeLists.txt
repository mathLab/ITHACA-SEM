SET(MESHUTILS_SOURCES
    CurveMesh.cpp
    Octant.cpp
    Octree.cpp
    SurfaceMesh.cpp
    SurfaceMeshing.cpp
    TriangleInterface.cpp
    )

SET(MESHUTILS_HEADERS
    CurvaturePoint.hpp
    MeshElem.hpp
    CurveMesh.h
    Octant.h
    Octree.h
    SurfaceMesh.h
    SurfaceMeshing.h
    TriangleInterface.h
    )

ADD_NEKTAR_LIBRARY(MeshUtils lib ${NEKTAR_LIBRARY_TYPE}
    ${MESHUTILS_SOURCES} ${MESHUTILS_HEADERS})
TARGET_LINK_LIBRARIES(MeshUtils LINK_PUBLIC LibUtilities)
TARGET_LINK_LIBRARIES(MeshUtils LINK_PRIVATE ${TRIANGLE_LIBRARY})

INSTALL(DIRECTORY ./
        DESTINATION ${NEKTAR_INCLUDE_DIR}/MeshUtils
        COMPONENT dev 
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
