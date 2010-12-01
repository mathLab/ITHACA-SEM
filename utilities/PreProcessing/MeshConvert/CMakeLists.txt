

SET(MeshConvertSources MeshConvert.cpp MeshElements.cpp ConvertGmsh.cpp
    ConvertPly.cpp Convert.cpp)
ADD_UTILITIES_EXECUTABLE(MeshConvert ${MeshConvertSources})

