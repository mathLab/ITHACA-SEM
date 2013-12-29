SET(FieldConvertHeaders 
    Module.h
    Field.hpp
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
ADD_NEKTAR_TEST(chan3D_tec_n10)
ADD_NEKTAR_TEST(chan3D_vtu)
ADD_NEKTAR_TEST(chan3D_vort)


IF (NEKTAR_USE_MPI)
    ADD_NEKTAR_TEST(chan3D_tec_par)
ENDIF (NEKTAR_USE_MPI)

