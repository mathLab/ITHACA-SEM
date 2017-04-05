SET(NekMeshHeaders
    InputModules/InputGmsh.h
    InputModules/InputNek.h
    InputModules/InputNek5000.h
    InputModules/InputNekpp.h
    InputModules/InputPly.h
    InputModules/InputSem.h
    InputModules/InputSwan.h
    InputModules/InputStarTec.h
    OutputModules/OutputGmsh.h
    OutputModules/OutputNekpp.h
    OutputModules/OutputSTL.h
    ProcessModules/ProcessBL.h
    ProcessModules/ProcessCurve.h
    ProcessModules/ProcessCurvedEdges.h
    ProcessModules/ProcessCyl.h
    ProcessModules/ProcessDetectSurf.h
    ProcessModules/ProcessExtractSurf.h
    ProcessModules/ProcessExtractTetPrismInterface.h
    ProcessModules/ProcessJac.h
    ProcessModules/ProcessLinkCheck.h
    ProcessModules/ProcessLinear.h
    ProcessModules/ProcessPerAlign.h
    ProcessModules/ProcessScalar.h
    ProcessModules/ProcessSpherigon.h
    ProcessModules/ProcessTetSplit.h
    ProcessModules/ProcessOptiExtract.h
    ProcessModules/ProcessInsertSurface.h
    ProcessModules/ProcessExtrude.h
    ProcessModules/ProcessVarOpti/ProcessVarOpti.h
    ProcessModules/ProcessVarOpti/NodeOpti.h
    ProcessModules/ProcessVarOpti/ElUtil.h
    ProcessModules/ProcessVarOpti/Evaluator.hxx
    ProcessModules/ProcessVarOpti/Hessian.hxx
)

SET(NekMeshSources
    InputModules/InputGmsh.cpp
    InputModules/InputNek.cpp
    InputModules/InputNek5000.cpp
    InputModules/InputNekpp.cpp
    InputModules/InputPly.cpp
    InputModules/InputSem.cpp
    InputModules/InputSwan.cpp
    InputModules/InputStarTec.cpp
    NekMesh.cpp
    OutputModules/OutputGmsh.cpp
    OutputModules/OutputNekpp.cpp
    OutputModules/OutputSTL.cpp
    ProcessModules/ProcessBL.cpp
    ProcessModules/ProcessCurve.cpp
    ProcessModules/ProcessCurvedEdges.cpp
    ProcessModules/ProcessCyl.cpp
    ProcessModules/ProcessDetectSurf.cpp
    ProcessModules/ProcessExtractSurf.cpp
    ProcessModules/ProcessExtractTetPrismInterface.cpp
    ProcessModules/ProcessJac.cpp
    ProcessModules/ProcessLinkCheck.cpp
    ProcessModules/ProcessLinear.cpp
    ProcessModules/ProcessPerAlign.cpp
    ProcessModules/ProcessScalar.cpp
    ProcessModules/ProcessSpherigon.cpp
    ProcessModules/ProcessTetSplit.cpp
    ProcessModules/ProcessOptiExtract.cpp
    ProcessModules/ProcessInsertSurface.cpp
    ProcessModules/ProcessExtrude.cpp
    ProcessModules/ProcessVarOpti/ProcessVarOpti.cpp
    ProcessModules/ProcessVarOpti/PreProcessing.cpp
    ProcessModules/ProcessVarOpti/NodeOpti.cpp
    ProcessModules/ProcessVarOpti/ElUtil.cpp
)

IF (NEKTAR_USE_CCM)
    SET(NekMeshHeaders ${NekMeshHeaders} InputModules/InputStar.h)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputStar.cpp)
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    SET(NekMeshHeaders ${NekMeshHeaders} InputModules/InputVtk.h OutputModules/OutputVtk.h)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputVtk.cpp OutputModules/OutputVtk.cpp)
ENDIF (NEKTAR_USE_VTK)

IF (NEKTAR_USE_MESHGEN)
    SET(NekMeshHeaders ${NekMeshHeaders}
            ProcessModules/ProcessVarOpti/NodeOptiCAD.h
            InputModules/InputMCF.h)
    SET(NekMeshSources ${NekMeshSources}
            ProcessModules/ProcessVarOpti/NodeOptiCAD.cpp
            InputModules/InputMCF.cpp)
ENDIF (NEKTAR_USE_MESHGEN)

# Don't use macro -- this way we can just link against NekMeshUtils and
# LocalRegions.
ADD_EXECUTABLE(NekMesh ${NekMeshSources})
SET_COMMON_PROPERTIES(NekMesh)
TARGET_LINK_LIBRARIES(NekMesh NekMeshUtils FieldUtils)
SET_PROPERTY(TARGET NekMesh PROPERTY FOLDER util)
INSTALL(TARGETS NekMesh
    RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT util)

IF (NEKTAR_USE_CCM)
   TARGET_LINK_LIBRARIES(NekMesh ${CCMIO_LIBRARIES})
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(NekMesh vtkCommon vtkIO vtkGraphics)
    ELSE ()
        TARGET_LINK_LIBRARIES(NekMesh vtkCommonCore vtkIOLegacy)
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)

# Nektar++
ADD_NEKTAR_TEST        (Nektar++/InvalidTetFace)
ADD_NEKTAR_TEST        (Nektar++/InvalidQuads)
ADD_NEKTAR_TEST        (Nektar++/Tube45Refinement)
ADD_NEKTAR_TEST        (Nektar++/Tube45Refinement_extractsurf)
ADD_NEKTAR_TEST        (Nektar++/CylQuadBl)
# Gmsh tests
ADD_NEKTAR_TEST        (Gmsh/CubeAllElements)
ADD_NEKTAR_TEST        (Gmsh/CubeHex)
ADD_NEKTAR_TEST        (Gmsh/CubeHexLinear)
ADD_NEKTAR_TEST        (Gmsh/CubePrism)
ADD_NEKTAR_TEST        (Gmsh/CubePrismLinear)
ADD_NEKTAR_TEST        (Gmsh/CubeTet)
ADD_NEKTAR_TEST        (Gmsh/CubeTetLinear)
IF (WIN32)
    ADD_NEKTAR_TEST    (Gmsh/Scalar_Windows)
ELSE ()
    ADD_NEKTAR_TEST    (Gmsh/Scalar)
ENDIF ()
ADD_NEKTAR_TEST        (Gmsh/SquareQuad)
ADD_NEKTAR_TEST        (Gmsh/SquareQuadLinear)
ADD_NEKTAR_TEST        (Gmsh/SquareTri)
ADD_NEKTAR_TEST        (Gmsh/SquareTriLinear)
# Nektar tests
ADD_NEKTAR_TEST        (Nektar/HexLinear)
ADD_NEKTAR_TEST        (Nektar/Tube45)
# Nek5000 tests
ADD_NEKTAR_TEST        (Nek5000/3da)
ADD_NEKTAR_TEST        (Nek5000/r1854a)
# StarCCM tests
ADD_NEKTAR_TEST        (StarTec/CubePer)
ADD_NEKTAR_TEST_LENGTHY(StarTec/StraightRW)

# MeshGen test
IF(NEKTAR_USE_MESHGEN)
    ADD_NEKTAR_TEST        (MeshGen/cylinder)
    ADD_NEKTAR_TEST        (MeshGen/sphere)
    ADD_NEKTAR_TEST        (MeshGen/2d-cad)
    ADD_NEKTAR_TEST        (MeshGen/2d-naca)
    ADD_NEKTAR_TEST        (MeshGen/t106step)
    ADD_NEKTAR_TEST_LENGTHY       (MeshGen/rev)
ENDIF()

# varopti tests
ADD_NEKTAR_TEST_LENGTHY       (MeshGen/varopti_cubesphere)
