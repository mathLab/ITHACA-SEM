SET(MeshConvertHeaders
    InputModules/InputGmsh.h
    InputModules/InputNek.h
    InputModules/InputNekpp.h
    InputModules/InputPly.h
    InputModules/InputSem.h
    InputModules/InputSwan.h
    InputModules/InputStarTec.h
    MeshElements.h
    Module.h
    OutputModules/OutputGmsh.h
    OutputModules/OutputNekpp.h
    ProcessModules/ProcessBL.h
    ProcessModules/ProcessCyl.h
    ProcessModules/ProcessDetectSurf.h
    ProcessModules/ProcessExtractSurf.h
    ProcessModules/ProcessExtractTetPrismInterface.h
    ProcessModules/ProcessJac.h
    ProcessModules/ProcessLinear.h
    ProcessModules/ProcessPerAlign.h
    ProcessModules/ProcessScalar.h
    ProcessModules/ProcessSpherigon.h
    ProcessModules/ProcessTetSplit.h
)

SET(MeshConvertSources
    InputModules/InputGmsh.cpp
    InputModules/InputNek.cpp
    InputModules/InputNekpp.cpp
    InputModules/InputPly.cpp
    InputModules/InputSem.cpp
    InputModules/InputSwan.cpp
    InputModules/InputStarTec.cpp
    MeshConvert.cpp
    MeshElements.cpp
    Module.cpp
    OutputModules/OutputGmsh.cpp
    OutputModules/OutputNekpp.cpp
    ProcessModules/ProcessBL.cpp
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
)

IF (NEKTAR_USE_CCM)
    SET(MeshConvertHeaders ${MeshConvertHeaders} InputModules/InputStar.h)
    SET(MeshConvertSources ${MeshConvertSources} InputModules/InputStar.cpp)
ENDIF (NEKTAR_USE_CCM)


IF (NEKTAR_USE_VTK)
    SET(MeshConvertHeaders ${MeshConvertHeaders} InputModules/InputVtk.h OutputModules/OutputVtk.h)
    SET(MeshConvertSources ${MeshConvertSources} InputModules/InputVtk.cpp OutputModules/OutputVtk.cpp)
ENDIF (NEKTAR_USE_VTK)

ADD_UTILITIES_EXECUTABLE(MeshConvert util ${MeshConvertSources} ${MeshConvertHeaders})

IF (NEKTAR_USE_CCM)
   TARGET_LINK_LIBRARIES(MeshConvert ccmio adf)
   SET_TARGET_PROPERTIES(MeshConvert PROPERTIES
    IMPORTED_LOCATION ${CCMIO_LIBRARY_PATH}/libccmio.a)
   INCLUDE_DIRECTORIES(MeshConvert ${CCMIO_INCLUDE_DIR})
   LINK_DIRECTORIES(${CCMIO_LIBRARY_PATH})
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(MeshConvert vtkCommon vtkIO vtkGraphics)
    ELSE ()
        TARGET_LINK_LIBRARIES(MeshConvert vtkCommonCore vtkIOLegacy)
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)

# Nektar++
ADD_NEKTAR_TEST        (Nektar++/InvalidTetFace)
# Gmsh tests
ADD_NEKTAR_TEST        (Gmsh/CubeAllElements)
ADD_NEKTAR_TEST        (Gmsh/CubeHex)
ADD_NEKTAR_TEST        (Gmsh/CubeHexLinearise)
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
