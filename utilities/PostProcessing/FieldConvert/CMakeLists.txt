SET(FieldConvertHeaders 
    Module.h
    Field.h
)

SET(FieldConvertSources 
    Module.cpp
    FieldConvert.cpp
)
  
ADD_UTILITIES_EXECUTABLE(FieldConvert util ${FieldConvertSources} ${FieldConvertHeaders})

