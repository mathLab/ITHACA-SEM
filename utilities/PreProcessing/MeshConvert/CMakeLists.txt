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
    ProcessDetectSurf.h
    ProcessExtractSurf.h
    ProcessJac.h
    ProcessPerAlign.h
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
    ProcessDetectSurf.cpp
    ProcessExtractSurf.cpp
    ProcessJac.cpp
    ProcessPerAlign.cpp
    ProcessSpherigon.cpp
    ProcessTetSplit.cpp
)

IF (NEKTAR_USE_VTK)
    SET(MeshConvertHeaders ${MeshConvertHeaders} InputVtk.h OutputVtk.h)
    SET(MeshConvertSources ${MeshConvertSources} InputVtk.cpp OutputVtk.cpp)
ENDIF (NEKTAR_USE_VTK)

ADD_UTILITIES_EXECUTABLE(MeshConvert util ${MeshConvertSources} ${MeshConvertHeaders})

#ADD_NEKTAR_TEST(MeshConvert_CubePer)
#ADD_NEKTAR_TEST_LENGTHY(MeshConvert_StraightRW)
