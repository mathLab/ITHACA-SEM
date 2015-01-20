SET(MeshConvertHeaders 
    InputGmsh.h
    InputNek.h
    InputNekpp.h
    InputPly.h
    InputSem.h
    InputSwan.h
    InputStarTec.h
    MeshElements.h
    Module.h
    OutputGmsh.h
    OutputNekpp.h
    ProcessBL.h
    ProcessCyl.h
    ProcessDetectSurf.h
    ProcessExtractSurf.h
    ProcessJac.h
    ProcessPerAlign.h
    ProcessScalar.h
    ProcessSpherigon.h
    ProcessTetSplit.h
)

SET(MeshConvertSources 
    InputGmsh.cpp
    InputNek.cpp
    InputNekpp.cpp
    InputPly.cpp
    InputSem.cpp
    InputSwan.cpp
    InputStarTec.cpp
    MeshConvert.cpp
    MeshElements.cpp
    Module.cpp
    OutputGmsh.cpp
    OutputNekpp.cpp
    ProcessBL.cpp
    ProcessCyl.cpp
    ProcessDetectSurf.cpp
    ProcessExtractSurf.cpp
    ProcessJac.cpp
    ProcessPerAlign.cpp
    ProcessScalar.cpp
    ProcessSpherigon.cpp
    ProcessTetSplit.cpp
)

IF (NEKTAR_USE_VTK)
    SET(MeshConvertHeaders ${MeshConvertHeaders} InputVtk.h OutputVtk.h)
    SET(MeshConvertSources ${MeshConvertSources} InputVtk.cpp OutputVtk.cpp)
ENDIF (NEKTAR_USE_VTK)

ADD_UTILITIES_EXECUTABLE(MeshConvert util ${MeshConvertSources} ${MeshConvertHeaders})

IF (NEKTAR_USE_VTK)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(MeshConvert vtkCommon vtkIO vtkGraphics)
    ELSE ()
        TARGET_LINK_LIBRARIES(MeshConvert vtkCommonCore vtkIOLegacy)
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)

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
