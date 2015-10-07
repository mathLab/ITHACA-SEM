SET(FieldConvertHeaders
    Module.h
    Field.hpp
    InputModules/InputDat.h
    InputModules/InputFld.h
    InputModules/InputXml.h
    InputModules/InputPts.h
    OutputModules/OutputInfo.h
    OutputModules/OutputTecplot.h
    OutputModules/OutputVtk.h
    OutputModules/OutputFld.h
    OutputModules/OutputXml.h
    ProcessModules/ProcessAddFld.h
    ProcessModules/ProcessBoundaryExtract.h
    ProcessModules/ProcessConcatenateFld.h
    ProcessModules/ProcessDeform.h
    ProcessModules/ProcessDisplacement.h
    ProcessModules/ProcessEquiSpacedOutput.h
    ProcessModules/ProcessGrad.h
    ProcessModules/ProcessInterpField.h
    ProcessModules/ProcessInterpPoints.h
    ProcessModules/ProcessInterpPointDataToFld.h
    ProcessModules/ProcessIsoContour.h
    ProcessModules/ProcessJacobianEnergy.h
    ProcessModules/ProcessPrintFldNorms.h
    ProcessModules/ProcessScaleInFld.h
    ProcessModules/ProcessSurfDistance.h
    ProcessModules/ProcessVorticity.h
    ProcessModules/ProcessScalGrad.h
    ProcessModules/ProcessMultiShear.h
    ProcessModules/ProcessWSS.h
    ProcessModules/ProcessC0Projection.h
    ProcessModules/ProcessQCriterion.h
)

SET(FieldConvertSources
    Module.cpp
    FieldConvert.cpp
    InputModules/InputDat.cpp
    InputModules/InputFld.cpp
    InputModules/InputXml.cpp
    InputModules/InputPts.cpp
    OutputModules/OutputInfo.cpp
    OutputModules/OutputTecplot.cpp
    OutputModules/OutputVtk.cpp
    OutputModules/OutputFld.cpp
    OutputModules/OutputXml.cpp
    ProcessModules/ProcessAddFld.cpp
    ProcessModules/ProcessBoundaryExtract.cpp
    ProcessModules/ProcessConcatenateFld.cpp
    ProcessModules/ProcessDeform.cpp
    ProcessModules/ProcessDisplacement.cpp
    ProcessModules/ProcessEquiSpacedOutput.cpp
    ProcessModules/ProcessGrad.cpp
    ProcessModules/ProcessInterpField.cpp
    ProcessModules/ProcessInterpPoints.cpp
    ProcessModules/ProcessInterpPointDataToFld.cpp
    ProcessModules/ProcessIsoContour.cpp
    ProcessModules/ProcessJacobianEnergy.cpp
    ProcessModules/ProcessPrintFldNorms.cpp
    ProcessModules/ProcessScaleInFld.cpp
    ProcessModules/ProcessVorticity.cpp
    ProcessModules/ProcessScalGrad.cpp
    ProcessModules/ProcessSurfDistance.cpp
    ProcessModules/ProcessMultiShear.cpp
    ProcessModules/ProcessWSS.cpp
    ProcessModules/ProcessC0Projection.cpp
    ProcessModules/ProcessQCriterion.cpp
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

