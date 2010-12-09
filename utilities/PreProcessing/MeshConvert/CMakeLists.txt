
SET(MeshConvertHeaders Convert.h ConvertGmsh.h ConvertPly.h MeshElements.h)

SET(MeshConvertSources MeshConvert.cpp MeshElements.cpp ConvertGmsh.cpp
    ConvertPly.cpp Convert.cpp)
ADD_UTILITIES_EXECUTABLE(MeshConvert ${MeshConvertSources} ${MeshConvertHeaders})

