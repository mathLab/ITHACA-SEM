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
)

IF (NEKTAR_USE_CCM)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputStar.cpp)
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputVtk.cpp OutputModules/OutputVtk.cpp)
ENDIF (NEKTAR_USE_VTK)

IF (NEKTAR_USE_MESHGEN)
    SET(NekMeshSources ${NekMeshSources}
            ProcessModules/ProcessVarOpti/NodeOptiCAD.cpp
            ProcessModules/ProcessVarOpti/ProcessVarOpti.cpp
            ProcessModules/ProcessVarOpti/PreProcessing.cpp
            ProcessModules/ProcessVarOpti/NodeOpti.cpp
            ProcessModules/ProcessVarOpti/ElUtil.cpp
            InputModules/InputCAD.cpp)
ENDIF (NEKTAR_USE_MESHGEN)

# Don't use macro -- this way we can just link against NekMeshUtils and
# LocalRegions.
ADD_UTILITIES_EXECUTABLE(NekMesh COMPONENT nekmesh
    SOURCES ${NekMeshSources}
    DEPENDS NekMeshUtils FieldUtils)

FINALISE_CPACK_COMPONENT(nekmesh DESCRIPTION "Nektar++ NekMesh utility")

IF (NEKTAR_USE_CCM)
   TARGET_LINK_LIBRARIES(NekMesh ${CCMIO_LIBRARIES})
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(NekMesh LINK_PUBLIC vtkCommon vtkIO vtkGraphics)
    ELSE ()
        TARGET_LINK_LIBRARIES(NekMesh LINK_PUBLIC vtkCommonCore vtkIOLegacy)
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)

# Nektar++
ADD_NEKTAR_TEST(Nektar++/InvalidTetFace)
ADD_NEKTAR_TEST(Nektar++/InvalidQuads)
ADD_NEKTAR_TEST(Nektar++/Tube45Refinement)
ADD_NEKTAR_TEST(Nektar++/Tube45Refinement_extractsurf)
ADD_NEKTAR_TEST(Nektar++/CylQuadBl)
# Gmsh tests
ADD_NEKTAR_TEST(Gmsh/CubeAllElements)
ADD_NEKTAR_TEST(Gmsh/CubeHex)
ADD_NEKTAR_TEST(Gmsh/CubeHexLinear)
ADD_NEKTAR_TEST(Gmsh/CubePrism)
ADD_NEKTAR_TEST(Gmsh/CubePrismLinear)
ADD_NEKTAR_TEST(Gmsh/CubeTet)
ADD_NEKTAR_TEST(Gmsh/CubeTetLinear)
ADD_NEKTAR_TEST(Gmsh/SquareQuad)
ADD_NEKTAR_TEST(Gmsh/SquareQuadLinear)
ADD_NEKTAR_TEST(Gmsh/SquareTri)
ADD_NEKTAR_TEST(Gmsh/SquareTriLinear)
# Nektar tests
ADD_NEKTAR_TEST(Nektar/HexLinear)
ADD_NEKTAR_TEST(Nektar/Tube45)
# Nek5000 tests
ADD_NEKTAR_TEST(Nek5000/3da)
ADD_NEKTAR_TEST(Nek5000/r1854a)
# StarCCM tests
ADD_NEKTAR_TEST(StarTec/CubePer)
ADD_NEKTAR_TEST(StarTec/StraightRW LENGTHY)

# MeshGen test
IF(NEKTAR_USE_MESHGEN)
    ADD_NEKTAR_TEST        (MeshGen/cylinder)
    ADD_NEKTAR_TEST        (MeshGen/sphere)
    ADD_NEKTAR_TEST        (MeshGen/2d-cad)
    ADD_NEKTAR_TEST        (MeshGen/2d-naca)
    # varopti tests
    ADD_NEKTAR_TEST        (MeshGen/varopti_cubesphere)
ENDIF()

# This has some issues due to hashing
IF (WIN32)
    ADD_NEKTAR_TEST(Gmsh/Scalar_Windows)
ELSE ()
    ADD_NEKTAR_TEST(Gmsh/Scalar)
ENDIF ()
