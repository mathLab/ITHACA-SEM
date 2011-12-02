
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

ADD_UTILITIES_EXECUTABLE(MeshConvert util ${MeshConvertSources} ${MeshConvertHeaders})

