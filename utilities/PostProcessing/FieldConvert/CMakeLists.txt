SET(FieldConvertHeaders 
    Module.h
    Field.h
    InputFld.h
    OutputTecplot.h
    OutputVtk.h
    OutputFld.h
    ProcessVorticity.h
)

SET(FieldConvertSources 
    Module.cpp
    FieldConvert.cpp
    InputFld.cpp
    OutputTecplot.cpp
    OutputVtk.cpp
    OutputFld.cpp
    ProcessVorticity.cpp
)
  
ADD_UTILITIES_EXECUTABLE(FieldConvert util ${FieldConvertSources} ${FieldConvertHeaders})


    ADD_NEKTAR_TEST(chan3D_tec)


