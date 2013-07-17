SET(FieldConvertHeaders 
    Module.h
    Field.h
    InputFld.h
    OutputVtk.h
    ProcessVorticity.h
)

SET(FieldConvertSources 
    Module.cpp
    FieldConvert.cpp
    InputFld.cpp
    OutputVtk.cpp
    ProcessVorticity.cpp
)
  
ADD_UTILITIES_EXECUTABLE(FieldConvert util ${FieldConvertSources} ${FieldConvertHeaders})

