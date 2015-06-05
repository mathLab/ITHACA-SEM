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
    OutputXml.h
    ProcessAddFld.h
    ProcessBoundaryExtract.h
    ProcessConcatenateFld.h
    ProcessDeform.h
    ProcessDisplacement.h
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
    OutputXml.cpp
    ProcessAddFld.cpp
    ProcessBoundaryExtract.cpp
    ProcessConcatenateFld.cpp
    ProcessDeform.cpp
    ProcessDisplacement.cpp
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
TARGET_LINK_LIBRARIES(FieldConvert SolverUtils)

#ADD_NEKTAR_TEST(chan3D_tec)
#ADD_NEKTAR_TEST(chan3D_tec_n10)
#ADD_NEKTAR_TEST(chan3D_vtu)
ADD_NEKTAR_TEST(chan3D_vort)
#ADD_NEKTAR_TEST(bfs_tec)
#ADD_NEKTAR_TEST(bfs_tec_rng)
ADD_NEKTAR_TEST(bfs_vort)
ADD_NEKTAR_TEST(bfs_vort_rng)
# ADD_NEKTAR_TEST(chan3D_pts)

# windows produces slightly differently formatted files which results in
# different hashes
# disabled these tests completely as they are system-dependent -- cc
#IF(WIN32)
#    ADD_NEKTAR_TEST(chan3D_probe_win)
#    ADD_NEKTAR_TEST(chan3D_equispacedoutput_win)
#    ADD_NEKTAR_TEST(chan3D_isocontour_win)
#ELSE(WIN32)
#    ADD_NEKTAR_TEST(chan3D_probe)
#    ADD_NEKTAR_TEST(chan3D_equispacedoutput)
#    ADD_NEKTAR_TEST(chan3D_isocontour)
#ENDIF(WIN32)


IF (NEKTAR_USE_MPI)
    #ADD_NEKTAR_TEST(chan3D_tec_par)
    #ADD_NEKTAR_TEST(chan3D_vtu_par)
    ADD_NEKTAR_TEST(chan3D_vort_par)
ENDIF (NEKTAR_USE_MPI)

