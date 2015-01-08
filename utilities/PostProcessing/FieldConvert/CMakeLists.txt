SET(FieldConvertHeaders 
    Module.h
    Field.hpp
    InputDat.h
    InputFld.h
    InputXml.h
    InputPts.h
    OutputInfo.h
    OutputTecplot.h
    OutputVtk.h
    OutputFld.h
    ProcessAddFld.h
    ProcessBoundaryExtract.h
    ProcessConcatenateFld.h
    ProcessEquiSpacedOutput.h
    ProcessInterpField.h
    ProcessInterpPoints.h
    ProcessInterpPointDataToFld.h
    ProcessIsoContour.h
    ProcessJacobianEnergy.h
    ProcessScaleInFld.h
    ProcessVorticity.h
    ProcessC0Projection.h
    ProcessQCriterion.h
)

SET(FieldConvertSources 
    Module.cpp
    FieldConvert.cpp
    InputDat.cpp
    InputFld.cpp
    InputXml.cpp
    InputPts.cpp
    OutputInfo.cpp
    OutputTecplot.cpp
    OutputVtk.cpp
    OutputFld.cpp
    ProcessAddFld.cpp
    ProcessBoundaryExtract.cpp
    ProcessConcatenateFld.cpp
    ProcessEquiSpacedOutput.cpp
    ProcessInterpField.cpp
    ProcessInterpPoints.cpp
    ProcessInterpPointDataToFld.cpp
    ProcessIsoContour.cpp
    ProcessJacobianEnergy.cpp
    ProcessScaleInFld.cpp
    ProcessVorticity.cpp
    ProcessC0Projection.cpp
    ProcessQCriterion.cpp
)
  
ADD_UTILITIES_EXECUTABLE(FieldConvert util ${FieldConvertSources} ${FieldConvertHeaders})


#ADD_NEKTAR_TEST(chan3D_tec)
#ADD_NEKTAR_TEST(chan3D_tec_n10)
#ADD_NEKTAR_TEST(chan3D_vtu)
ADD_NEKTAR_TEST(chan3D_vort)
#ADD_NEKTAR_TEST(bfs_tec)
#ADD_NEKTAR_TEST(bfs_tec_rng)
ADD_NEKTAR_TEST(bfs_vort)
ADD_NEKTAR_TEST(bfs_vort_rng)


IF (NEKTAR_USE_MPI)
    #ADD_NEKTAR_TEST(chan3D_tec_par)
    #ADD_NEKTAR_TEST(chan3D_vtu_par)
    ADD_NEKTAR_TEST(chan3D_vort_par)
ENDIF (NEKTAR_USE_MPI)

