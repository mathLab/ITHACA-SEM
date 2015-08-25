SET(MESHUTILS_SOURCES
    SurfaceMeshing/CurveMesh.cpp
    Octree/Octant.cpp
    Octree/Octree.cpp
    SurfaceMeshing/SurfaceMesh.cpp
    SurfaceMeshing/SurfaceMeshing.cpp
    SurfaceMeshing/OptimiseFunctions.cpp
    SurfaceMeshing/OptimiseUtils.cpp
    ExtLibInterface/TriangleInterface.cpp
    ExtLibInterface/TetGenInterface.cpp
    TetMeshing/TetMesh.cpp
    )

SET(MESHUTILS_HEADERS
    Octree/CurvaturePoint.hpp
    MeshElem.hpp
    SurfaceMeshing/CurveMesh.h
    Octree/Octant.h
    Octree/Octree.h
    SurfaceMeshing/SurfaceMesh.h
    SurfaceMeshing/SurfaceMeshing.h
    ExtLibInterface/TriangleInterface.h
    ExtLibInterface/TetGenInterface.h
    TetMeshing/TetMesh.h
    )

ADD_NEKTAR_LIBRARY(MeshUtils lib ${NEKTAR_LIBRARY_TYPE}
    ${MESHUTILS_SOURCES} ${MESHUTILS_HEADERS})
TARGET_LINK_LIBRARIES(MeshUtils LINK_PUBLIC LibUtilities)
TARGET_LINK_LIBRARIES(MeshUtils LINK_PRIVATE ${TRIANGLE_LIBRARY})
TARGET_LINK_LIBRARIES(MeshUtils LINK_PRIVATE ${TETGEN_LIBRARY})

INSTALL(DIRECTORY ./
        DESTINATION ${NEKTAR_INCLUDE_DIR}/MeshUtils
        COMPONENT dev
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
