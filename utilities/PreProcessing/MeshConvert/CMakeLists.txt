
SET(MeshConvertHeaders 
    Module.h
    InputGmsh.h
    InputNek.h
    InputPly.h
    InputSem.h
    OutputNekpp.h
    MeshElements.h
)

SET(MeshConvertSources 
    Module.cpp
    InputGmsh.cpp
    InputNek.cpp
    InputPly.cpp
    InputSem.cpp
    OutputNekpp.cpp
    MeshConvert.cpp
    MeshElements.cpp)

IF (NEKTAR_USE_VTK)
    SET(MeshConvertHeaders ${MeshConvertHeaders} InputVtk.h OutputVtk.h)
    SET(MeshConvertSources ${MeshConvertSources} InputVtk.cpp OutputVtk.cpp)
ENDIF (NEKTAR_USE_VTK)

ADD_UTILITIES_EXECUTABLE(MeshConvert util ${MeshConvertSources} ${MeshConvertHeaders})

