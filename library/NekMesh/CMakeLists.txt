SET(NEKMESH_SOURCES
    Module/Module.cpp
    Module/InputModules/InputGmsh.cpp
    Module/InputModules/InputNek.cpp
    Module/InputModules/InputNek5000.cpp
    Module/InputModules/InputNekpp.cpp
    Module/InputModules/InputPly.cpp
    Module/InputModules/InputSem.cpp
    Module/InputModules/InputSwan.cpp
    Module/InputModules/InputStarTec.cpp
    Module/OutputModules/OutputGmsh.cpp
    Module/OutputModules/OutputNekpp.cpp
    Module/OutputModules/OutputSTL.cpp
    Module/ProcessModules/ProcessBL.cpp
    Module/ProcessModules/ProcessCurve.cpp
    Module/ProcessModules/ProcessCurvedEdges.cpp
    Module/ProcessModules/ProcessCyl.cpp
    Module/ProcessModules/ProcessDetectSurf.cpp
    Module/ProcessModules/ProcessExtractSurf.cpp
    Module/ProcessModules/ProcessExtractTetPrismInterface.cpp
    Module/ProcessModules/ProcessJac.cpp
    Module/ProcessModules/ProcessLinkCheck.cpp
    Module/ProcessModules/ProcessLinear.cpp
    Module/ProcessModules/ProcessPerAlign.cpp
    Module/ProcessModules/ProcessScalar.cpp
    Module/ProcessModules/ProcessSpherigon.cpp
    Module/ProcessModules/ProcessTetSplit.cpp
    Module/ProcessModules/ProcessOptiExtract.cpp
    Module/ProcessModules/ProcessInsertSurface.cpp
    Module/ProcessModules/ProcessExtrude.cpp
    Module/ProcessModules/ProcessVarOpti/ProcessVarOpti.cpp
    Module/ProcessModules/ProcessVarOpti/PreProcessing.cpp
    Module/ProcessModules/ProcessVarOpti/NodeOpti.cpp
    Module/ProcessModules/ProcessVarOpti/ElUtil.cpp

    MeshElements/Mesh.cpp
    MeshElements/Element.cpp
    MeshElements/Edge.cpp
    MeshElements/Face.cpp
    MeshElements/BooleanOperators.cpp
    MeshElements/Point.cpp
    MeshElements/Node.cpp
    MeshElements/Line.cpp
    MeshElements/Triangle.cpp
    MeshElements/Quadrilateral.cpp
    MeshElements/Tetrahedron.cpp
    MeshElements/Pyramid.cpp
    MeshElements/Prism.cpp
    MeshElements/Hexahedron.cpp
    CADSystem/CADSystem.cpp
    CADSystem/CADSurf.cpp
    CADSystem/CADVert.cpp
    CADSystem/CADCurve.cpp
    CADSystem/ProcessLoadCAD.cpp
    )

SET(NEKMESH_HEADERS
    NekMeshDeclspec.h
    Module/Module.h
    Module/InputModules/InputGmsh.h
    Module/InputModules/InputNek.h
    Module/InputModules/InputNek5000.h
    Module/InputModules/InputNekpp.h
    Module/InputModules/InputPly.h
    Module/InputModules/InputSem.h
    Module/InputModules/InputSwan.h
    Module/InputModules/InputStarTec.h
    Module/OutputModules/OutputGmsh.h
    Module/OutputModules/OutputNekpp.h
    Module/OutputModules/OutputSTL.h
    Module/ProcessModules/ProcessBL.h
    Module/ProcessModules/ProcessCurve.h
    Module/ProcessModules/ProcessCurvedEdges.h
    Module/ProcessModules/ProcessCyl.h
    Module/ProcessModules/ProcessDetectSurf.h
    Module/ProcessModules/ProcessExtractSurf.h
    Module/ProcessModules/ProcessExtractTetPrismInterface.h
    Module/ProcessModules/ProcessJac.h
    Module/ProcessModules/ProcessLinkCheck.h
    Module/ProcessModules/ProcessLinear.h
    Module/ProcessModules/ProcessPerAlign.h
    Module/ProcessModules/ProcessScalar.h
    Module/ProcessModules/ProcessSpherigon.h
    Module/ProcessModules/ProcessTetSplit.h
    Module/ProcessModules/ProcessOptiExtract.h
    Module/ProcessModules/ProcessInsertSurface.h
    Module/ProcessModules/ProcessExtrude.h
    Module/ProcessModules/ProcessVarOpti/ProcessVarOpti.h
    Module/ProcessModules/ProcessVarOpti/NodeOpti.h
    Module/ProcessModules/ProcessVarOpti/ElUtil.h
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
    MeshElements/HOAlignment.h
    MeshElements/ElementConfig.h
    CADSystem/CADObject.h
    CADSystem/CADSystem.h
    CADSystem/CADVert.h
    CADSystem/CADCurve.h
    CADSystem/CADSurf.h
    CADSystem/ProcessLoadCAD.h
    )

IF(NEKTAR_USE_MESHGEN)
    SET(NEKMESH_SOURCES ${NEKMESH_SOURCES}
        SurfaceMeshing/CurveMesh.cpp
        Octree/Octant.cpp
        Octree/Octree.cpp
        Octree/ProcessLoadOctree.cpp
        SurfaceMeshing/SurfaceMesh.cpp
        SurfaceMeshing/HOSurfaceMesh.cpp
        SurfaceMeshing/FaceMesh.cpp
        SurfaceMeshing/OptimiseFunctions.cpp
        ExtLibInterface/TriangleInterface.cpp
        ExtLibInterface/TetGenInterface.cpp
        VolumeMeshing/VolumeMesh.cpp
        VolumeMeshing/TetMeshing/TetMesh.cpp
        VolumeMeshing/BLMeshing/BLMesh.cpp
        Optimisation/BGFS-B.cpp
        CADSystem/OCE/CADSystemOCE.cpp
        CADSystem/OCE/CADVertOCE.cpp
        CADSystem/OCE/CADCurveOCE.cpp
        CADSystem/OCE/CADSurfOCE.cpp
        CADSystem/OCE/TransfiniteSurface.cpp
        2DGenerator/2DGenerator.cpp
        Module/ProcessModules/ProcessVarOpti/NodeOptiCAD.cpp
        Module/ProcessModules/ProcessProjectCAD.cpp
        Module/InputModules/InputMCF.cpp
    )

    SET(NEKMESH_HEADERS ${NEKMESH_HEADERS}
        Octree/SourcePoint.hpp
        SurfaceMeshing/CurveMesh.h
        Octree/Octant.h
        Octree/Octree.h
        Octree/ProcessLoadOctree.h
        SurfaceMeshing/SurfaceMesh.h
        SurfaceMeshing/HOSurfaceMesh.h
        SurfaceMeshing/FaceMesh.h
        SurfaceMeshing/OptimiseFunctions.h
        ExtLibInterface/TriangleInterface.h
        ExtLibInterface/TetGenInterface.h
        VolumeMeshing/VolumeMesh.h
        VolumeMeshing/TetMeshing/TetMesh.h
        VolumeMeshing/BLMeshing/BLMesh.h
        CADSystem/OCE/OpenCascade.h
        CADSystem/OCE/CADSystemOCE.h
        CADSystem/OCE/CADVertOCE.h
        CADSystem/OCE/CADCurveOCE.h
        CADSystem/OCE/CADSurfOCE.h
        CADSystem/OCE/GeoParser.hpp
        CADSystem/OCE/TransfiniteSurface.h
        Optimisation/BGFS-B.h
        Optimisation/OptimiseObj.h
        2DGenerator/2DGenerator.h
        Module/ProcessModules/ProcessVarOpti/NodeOptiCAD.h
        Module/ProcessModules/ProcessProjectCAD.h
        Module/InputModules/InputMCF.h
    )
ENDIF()

IF(NEKTAR_USE_CFI)
    ADD_DEFINITIONS(-DNEKTAR_USE_CFI)
    SET(NEKMESH_HEADERS ${NEKMESH_HEADERS}
        CADSystem/CFI/CADSystemCFI.h
        CADSystem/CFI/CADVertCFI.h
        CADSystem/CFI/CADCurveCFI.h
        CADSystem/CFI/CADSurfCFI.h
        CADSystem/CFI/CADElementCFI.h
        Module/InputModules/InputCADfix.h
        Module/OutputModules/OutputCADfix.h
    )
    SET(NEKMESH_SOURCES ${NEKMESH_SOURCES}
        CADSystem/CFI/CADSystemCFI.cpp
        CADSystem/CFI/CADVertCFI.cpp
        CADSystem/CFI/CADCurveCFI.cpp
        CADSystem/CFI/CADSurfCFI.cpp
        Module/InputModules/InputCADfix.cpp
        Module/OutputModules/OutputCADfix.cpp
    )
ENDIF()

IF (NEKTAR_USE_CCM)
    SET(NEKMESH_SOURCES ${NEKMESH_SOURCES} Module/InputModules/InputStar.cpp)
    SET(NEKMESH_HEADERS ${NEKMESH_HEADERS} Module/InputModules/InputStar.h)
ENDIF()

IF (NEKTAR_USE_VTK)
    SET(NEKMESH_SOURCES ${NEKMESH_SOURCES}
        Module/InputModules/InputVtk.cpp
        Module/OutputModules/OutputVtk.cpp)
    SET(NEKMESH_HEADERS ${NEKMESH_HEADERS}
        Module/InputModules/InputVtk.h
        Module/OutputModules/OutputVtk.h)
ENDIF()

ADD_DEFINITIONS(-DNEKMESH_EXPORTS)

ADD_NEKTAR_LIBRARY(NekMesh
    SOURCES ${NEKMESH_SOURCES}
    HEADERS ${NEKMESH_HEADERS}
    DEPENDS LocalRegions
    SUMMARY "Nektar++ NekMesh library"
    DESCRIPTION "This library provides mesh generation and routines for the NekMesh executable.")

IF(NEKTAR_USE_MESHGEN)
    TARGET_LINK_LIBRARIES(NekMesh LINK_PRIVATE ${TETGEN_LIBRARY})
    TARGET_LINK_LIBRARIES(NekMesh LINK_PRIVATE ${TRIANGLE_LIBRARY})
    TARGET_LINK_LIBRARIES(NekMesh LINK_PRIVATE ${OCC_LIBRARIES})

    ADD_DEPENDENCIES(NekMesh oce-0.17 tetgen-1.5 triangle-1.6)

    IF(NEKTAR_USE_CFI)
	TARGET_LINK_LIBRARIES(NekMesh LINK_PUBLIC ${CFI_LIBRARY_API})
	TARGET_LINK_LIBRARIES(NekMesh LINK_PUBLIC ${CFI_LIBRARY_CXX})
    ENDIF()
ENDIF()

IF (NEKTAR_USE_CCM)
   TARGET_LINK_LIBRARIES(NekMesh LINK_PRIVATE ${CCMIO_LIBRARIES})
   ADD_DEPENDENCIES(NekMesh libccmio-2.6.1)
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(NekMesh LINK_PUBLIC vtkCommon vtkIO vtkGraphics)
    ELSE ()
        TARGET_LINK_LIBRARIES(NekMesh LINK_PUBLIC vtkCommonCore vtkIOLegacy)
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)

IF (NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()

INSTALL(DIRECTORY ./
        DESTINATION ${NEKTAR_INCLUDE_DIR}/NekMesh
        COMPONENT dev
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
