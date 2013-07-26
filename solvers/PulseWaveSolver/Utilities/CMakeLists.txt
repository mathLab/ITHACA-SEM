SET(Fld2TecplotSources  Fld2Tecplot.cpp 
       ../EquationSystems/PulseWaveSystem.h
       ../EquationSystems/PulseWaveSystem.cpp
       ../EquationSystems/PulseWavePropagation.h
       ../EquationSystems/PulseWavePropagation.cpp
       )


ADD_SOLVER_EXECUTABLE(Fld2Tecplot solvers ${Fld2TecplotSources})

