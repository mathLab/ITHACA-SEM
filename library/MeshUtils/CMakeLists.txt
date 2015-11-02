SET(MESHUTILS_SOURCES
    MeshElements/MeshElements.cpp
    )
IF(NEKTAR_USE_MESH)
    SET(MESHUTILS_SOURCES ${MESHUTILS_SOURCES}
        SurfaceMeshing/CurveMesh.cpp
        Octree/Octant.cpp
        Octree/Octree.cpp
        SurfaceMeshing/SurfaceMesh.cpp
        SurfaceMeshing/SurfaceMeshHOMesh.cpp
        SurfaceMeshing/FaceMesh.cpp
        SurfaceMeshing/OptimiseFunctions.cpp
        SurfaceMeshing/OptimiseUtils.cpp
        ExtLibInterface/TriangleInterface.cpp
        ExtLibInterface/TetGenInterface.cpp
        TetMeshing/TetMesh.cpp
    )
ENDIF()

SET(MESHUTILS_HEADERS
    MeshElements/MeshElements.h
    MeshElements/Node.h
    MeshElements/Edge.h
    MeshElements/Face.h
    MeshElements/Element.h
    MeshElements/Composite.h
    MeshElements/Mesh.h
    MeshElements/Point.h
    MeshElements/Line.h
    MeshElements/Triangle.h
    MeshElements/Quadrilateral.h
    MeshElements/Tetrahedron.h
    MeshElements/Pyramid.h
    MeshElements/Prism.h
    MeshElements/Hexahedron.h
    )
IF(NEKTAR_USE_MESH)
    SET(MESHUTILS_HEADERS ${MESHUTILS_HEADERS}
        Octree/CurvaturePoint.hpp
        MeshElem.hpp
        SurfaceMeshing/CurveMesh.h
        Octree/Octant.h
        Octree/Octree.h
        SurfaceMeshing/SurfaceMesh.h
        SurfaceMeshing/FaceMesh.h
        ExtLibInterface/TriangleInterface.h
        ExtLibInterface/TetGenInterface.h
        TetMeshing/TetMesh.h
    )
ENDIF()

IF(NEKTAR_USE_OCC)
    SET(MESHUTILS_SOURCES ${MESHUTILS_SOURCES}
        CADSystem/CADSystem.cpp
        CADSystem/CADCurve.cpp
        CADSystem/CADSurf.cpp
        )

    SET(MESHUTILS_HEADERS ${MESHUTILS_HEADERS}
        CADSystem/CADSystem.h
        CADSystem/CADVert.h
        CADSystem/CADCurve.h
        CADSystem/CADSurf.h
        CADSystem/OpenCascade.h
        )

ENDIF(NEKTAR_USE_OCC)

ADD_NEKTAR_LIBRARY(MeshUtils lib ${NEKTAR_LIBRARY_TYPE}
    ${MESHUTILS_SOURCES} ${MESHUTILS_HEADERS})
TARGET_LINK_LIBRARIES(MeshUtils LINK_PUBLIC SolverUtils)

IF(NEKTAR_USE_MESH)
    TARGET_LINK_LIBRARIES(MeshUtils LINK_PRIVATE ${TRIANGLE_LIBRARY})
    TARGET_LINK_LIBRARIES(MeshUtils LINK_PRIVATE ${TETGEN_LIBRARY})
ENDIF()

IF(NEKTAR_USE_OCC)
    SET(OCC_DEF LIN LININTEL HAVE_WOK_CONFIG_H HAVE_CONFIG_H CSFDB)
    SET_TARGET_PROPERTIES(MeshUtils PROPERTIES COMPILE_DEFINITIONS "${OCC_DEF}")
    TARGET_LINK_LIBRARIES(MeshUtils LINK_PRIVATE ${OCC_LIBS})
    ADD_DEPENDENCIES(MeshUtils opencascade-6.8)
ENDIF(NEKTAR_USE_OCC)

INSTALL(DIRECTORY ./
        DESTINATION ${NEKTAR_INCLUDE_DIR}/MeshUtils
        COMPONENT dev
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
