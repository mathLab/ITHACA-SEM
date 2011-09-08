
SET(MeshConvertHeaders 
    Convert.h 
    ConvertGmsh.h 
    ConvertPly.h 
    ConvertSem.h
    MeshElements.h
)

SET(MeshConvertSources 
    MeshConvert.cpp 
    MeshElements.cpp 
    ConvertGmsh.cpp
    ConvertPly.cpp 
    ConvertSem.cpp
    Convert.cpp)

ADD_UTILITIES_EXECUTABLE(MeshConvert util ${MeshConvertSources} ${MeshConvertHeaders})

