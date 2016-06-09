SET(NekMeshSources
    InputModules/InputGmsh.cpp
    InputModules/InputNek.cpp
    InputModules/InputNekpp.cpp
    InputModules/InputPly.cpp
    InputModules/InputSem.cpp
    InputModules/InputSwan.cpp
    InputModules/InputStarTec.cpp
    NekMesh.cpp
    Module.cpp
    OutputModules/OutputGmsh.cpp
    OutputModules/OutputNekpp.cpp
    ProcessModules/ProcessBL.cpp
    ProcessModules/ProcessCurve.cpp
    ProcessModules/ProcessCurvedEdges.cpp
    ProcessModules/ProcessCyl.cpp
    ProcessModules/ProcessDetectSurf.cpp
    ProcessModules/ProcessExtractSurf.cpp
    ProcessModules/ProcessExtractTetPrismInterface.cpp
    ProcessModules/ProcessJac.cpp
    ProcessModules/ProcessLinear.cpp
    ProcessModules/ProcessPerAlign.cpp
    ProcessModules/ProcessScalar.cpp
    ProcessModules/ProcessSpherigon.cpp
    ProcessModules/ProcessTetSplit.cpp
    ProcessModules/ProcessOptiExtract.cpp
)

IF (NEKTAR_USE_CCM)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputStar.cpp)
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputVtk.cpp OutputModules/OutputVtk.cpp)
ENDIF (NEKTAR_USE_VTK)

IF (NEKTAR_USE_MESHGEN)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputCAD.cpp)
ENDIF (NEKTAR_USE_MESHGEN)

# Don't use macro -- this way we can just link against NekMeshUtils and
# LocalRegions.
ADD_UTILITIES_EXECUTABLE(NekMesh COMPONENT nekmesh
    SOURCES ${NekMeshSources}
    DEPENDS NekMeshUtils LocalRegions)

IF (NEKTAR_USE_CCM)
   TARGET_LINK_LIBRARIES(NekMesh ccmio adf)
   SET_TARGET_PROPERTIES(NekMesh PROPERTIES
       IMPORTED_LOCATION ${CCMIO_LIBRARY_PATH}/libccmio.a)
   INCLUDE_DIRECTORIES(NekMesh ${CCMIO_INCLUDE_DIR})
   LINK_DIRECTORIES(${CCMIO_LIBRARY_PATH})
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
ADD_NEKTAR_TEST        (Nektar++/Tube45Refinement)
ADD_NEKTAR_TEST        (Nektar++/Tube45Refinement_extractsurf)
ADD_NEKTAR_TEST        (Nektar++/CylQuadBl)
# Gmsh tests
ADD_NEKTAR_TEST        (Gmsh/CubeAllElements)
ADD_NEKTAR_TEST        (Gmsh/CubeHex)
ADD_NEKTAR_TEST        (Gmsh/CubePrism)
ADD_NEKTAR_TEST        (Gmsh/CubeTet)
IF (WIN32)
    ADD_NEKTAR_TEST    (Gmsh/Scalar_Windows)
ELSE ()
    ADD_NEKTAR_TEST    (Gmsh/Scalar)
ENDIF ()
ADD_NEKTAR_TEST        (Gmsh/SquareQuad)
ADD_NEKTAR_TEST        (Gmsh/SquareTri)
# Nektar tests
ADD_NEKTAR_TEST        (Nektar/HexLinear)
ADD_NEKTAR_TEST        (Nektar/Tube45)
# StarCCM tests
ADD_NEKTAR_TEST        (StarTec/CubePer)
ADD_NEKTAR_TEST_LENGTHY(StarTec/StraightRW)
