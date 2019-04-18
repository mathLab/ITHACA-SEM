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
    SET(NekMeshSources ${NekMeshSources} InputModules/InputStar.cpp)
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputVtk.cpp OutputModules/OutputVtk.cpp)
ENDIF (NEKTAR_USE_VTK)

IF (NEKTAR_USE_MESHGEN)
    SET(NekMeshSources ${NekMeshSources}
            ProcessModules/ProcessVarOpti/NodeOptiCAD.cpp
            ProcessModules/ProcessProjectCAD.cpp
            InputModules/InputMCF.cpp)
ENDIF (NEKTAR_USE_MESHGEN)

# Don't use macro -- this way we can just link against NekMeshUtils and
# LocalRegions.
ADD_UTILITIES_EXECUTABLE(NekMesh COMPONENT nekmesh
    SOURCES ${NekMeshSources}
    DEPENDS NekMeshUtils LocalRegions)

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

# Nektar++
ADD_NEKTAR_TEST(Nektar++/InvalidTetFace)
ADD_NEKTAR_TEST(Nektar++/InvalidQuads)
ADD_NEKTAR_TEST(Nektar++/Tube45Refinement)
ADD_NEKTAR_TEST(Nektar++/Tube45Refinement_extractsurf)
ADD_NEKTAR_TEST(Nektar++/CylQuadBl)
ADD_NEKTAR_TEST(Nektar++/SurfaceExtractionCurvedEdge)
ADD_NEKTAR_TEST(Nektar++/SurfaceExtractionCurvedFace)
# Gmsh tests
ADD_NEKTAR_TEST(Gmsh/CubeAllElements_2_2)
ADD_NEKTAR_TEST(Gmsh/CubeHex_2_2)
ADD_NEKTAR_TEST(Gmsh/CubeHexLinear_2_2)
ADD_NEKTAR_TEST(Gmsh/CubePrism_2_2)
ADD_NEKTAR_TEST(Gmsh/CubePrismLinear_2_2)
ADD_NEKTAR_TEST(Gmsh/CubeTet_2_2)
ADD_NEKTAR_TEST(Gmsh/CubeTetLinear_2_2)
ADD_NEKTAR_TEST(Gmsh/SquareQuad_2_2)
ADD_NEKTAR_TEST(Gmsh/SquareQuadLinear_2_2)
ADD_NEKTAR_TEST(Gmsh/SquareTri_2_2)
ADD_NEKTAR_TEST(Gmsh/SquareTriLinear_2_2)
ADD_NEKTAR_TEST(Gmsh/CubeAllElements_4_0)
ADD_NEKTAR_TEST(Gmsh/CubeHex_4_0)
ADD_NEKTAR_TEST(Gmsh/CubeHexLinear_4_0)
ADD_NEKTAR_TEST(Gmsh/CubePrism_4_0)
ADD_NEKTAR_TEST(Gmsh/CubePrismLinear_4_0)
ADD_NEKTAR_TEST(Gmsh/CubeTet_4_0)
ADD_NEKTAR_TEST(Gmsh/CubeTetLinear_4_0)
ADD_NEKTAR_TEST(Gmsh/SquareCirc_4_0)
ADD_NEKTAR_TEST(Gmsh/SquareQuad_4_0)
ADD_NEKTAR_TEST(Gmsh/SquareQuadLinear_4_0)
ADD_NEKTAR_TEST(Gmsh/SquareTri_4_0)
ADD_NEKTAR_TEST(Gmsh/SquareTriLinear_4_0)
ADD_NEKTAR_TEST(Gmsh/CylinderRotPer2D)
# Nektar tests
ADD_NEKTAR_TEST(Nektar/HexLinear)
ADD_NEKTAR_TEST(Nektar/Tube45)
ADD_NEKTAR_TEST(Nektar/UKMesh)
# Nek5000 tests
ADD_NEKTAR_TEST(Nek5000/3da)
ADD_NEKTAR_TEST(Nek5000/r1854a)
# StarCCM tests
ADD_NEKTAR_TEST(StarTec/CubePer)
ADD_NEKTAR_TEST(StarTec/StraightRW LENGTHY)

# MeshGen test
IF(NEKTAR_USE_MESHGEN)
    ADD_NEKTAR_TEST(MeshGen/cylinder)
    ADD_NEKTAR_TEST(MeshGen/sphere)
    ADD_NEKTAR_TEST(MeshGen/2d-cad)
    ADD_NEKTAR_TEST(MeshGen/2d-naca)
    ADD_NEKTAR_TEST(MeshGen/geo)
    ADD_NEKTAR_TEST(MeshGen/t106c)
    ADD_NEKTAR_TEST(MeshGen/t106step)
    ADD_NEKTAR_TEST(MeshGen/rev LENGTHY)

    IF(NEKTAR_USE_CCM)              
        ADD_NEKTAR_TEST(StarCCM/Ahmed_body)
    ENDIF() 
ENDIF()

IF (NEKTAR_USE_HDF5)
    ADD_NEKTAR_TEST(Gmsh/CubeTet_4_0_hdf5)
ENDIF()

# varopti tests
ADD_NEKTAR_TEST(MeshGen/varopti_cubesphere LENGTHY)
