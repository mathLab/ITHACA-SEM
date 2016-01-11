SET(NEKMESHUTILS_SOURCES
    MeshElements/Mesh.cpp
    MeshElements/Element.cpp
    MeshElements/BooleanOperators.cpp
    MeshElements/Point.cpp
    MeshElements/Line.cpp
    MeshElements/Triangle.cpp
    MeshElements/Quadrilateral.cpp
    MeshElements/Tetrahedron.cpp
    MeshElements/Pyramid.cpp
    MeshElements/Prism.cpp
    MeshElements/Hexahedron.cpp
    )
IF(NEKTAR_USE_MESH)
    SET(NEKMESHUTILS_SOURCES ${NEKMESHUTILS_SOURCES}
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
        BLMeshing/BLMesh.cpp
        CADSystem/CADSystem.cpp
        CADSystem/CADCurve.cpp
        CADSystem/CADSurf.cpp
    )
ENDIF()

SET(NEKMESHUTILS_HEADERS
    NekMeshUtilsDeclspec.h
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
    SET(NEKMESHUTILS_HEADERS ${NEKMESHUTILS_HEADERS}
        Octree/CurvaturePoint.hpp
        SurfaceMeshing/CurveMesh.h
        Octree/Octant.h
        Octree/Octree.h
        SurfaceMeshing/SurfaceMesh.h
        SurfaceMeshing/FaceMesh.h
        ExtLibInterface/TriangleInterface.h
        ExtLibInterface/TetGenInterface.h
        TetMeshing/TetMesh.h
        BLMeshing/BLMesh.h
        CADSystem/CADSystem.h
        CADSystem/CADVert.h
        CADSystem/CADCurve.h
        CADSystem/CADSurf.h
        CADSystem/OpenCascade.h
    )
ENDIF()

ADD_NEKTAR_LIBRARY(NekMeshUtils lib ${NEKTAR_LIBRARY_TYPE}
    ${NEKMESHUTILS_SOURCES} ${NEKMESHUTILS_HEADERS})
TARGET_LINK_LIBRARIES(NekMeshUtils LINK_PUBLIC SolverUtils)

IF(NEKTAR_USE_MESH)
    TARGET_LINK_LIBRARIES(NekMeshUtils LINK_PRIVATE ${TRIANGLE_LIBRARY})
    TARGET_LINK_LIBRARIES(NekMeshUtils LINK_PRIVATE ${TETGEN_LIBRARY})
ENDIF()

IF(NEKTAR_USE_MESH)
    SET(OCC_DEF LIN LININTEL HAVE_WOK_CONFIG_H HAVE_CONFIG_H CSFDB)
    SET_TARGET_PROPERTIES(NekMeshUtils PROPERTIES COMPILE_DEFINITIONS "${OCC_DEF}")
    TARGET_LINK_LIBRARIES(NekMeshUtils LINK_PUBLIC ${OCC_LIBS})
    ADD_DEPENDENCIES(NekMeshUtils opencascade-6.8)
ENDIF(NEKTAR_USE_MESH)

IF (NEKTAR_USE_MESH)
    SET_TARGET_PROPERTIES(NekMeshUtils PROPERTIES COMPILE_DEFINITIONS "MESHGEN")
ENDIF()

INSTALL(DIRECTORY ./
        DESTINATION ${NEKTAR_INCLUDE_DIR}/NekMeshUtils
        COMPONENT dev
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
