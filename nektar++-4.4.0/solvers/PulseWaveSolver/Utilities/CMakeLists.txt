SET(Fld2TecplotSources  Fld2Tecplot.cpp 
       ../EquationSystems/PulseWaveSystem.h
       ../EquationSystems/PulseWaveSystem.cpp
       ../EquationSystems/PulseWavePropagation.h
       ../EquationSystems/PulseWavePropagation.cpp
       ../EquationSystems/PulseWaveBoundary.h
       ../EquationSystems/PulseWaveBoundary.cpp
       ../EquationSystems/PulseWavePressureArea.cpp
       ../EquationSystems/TimeDependentInflow.h
       ../EquationSystems/TimeDependentInflow.cpp
       ../EquationSystems/QInflow.h
       ../EquationSystems/QInflow.cpp
       ../EquationSystems/TerminalOutflow.h
       ../EquationSystems/TerminalOutflow.cpp
       ../EquationSystems/ROutflow.h
       ../EquationSystems/ROutflow.cpp
       ../EquationSystems/RCROutflow.h
       ../EquationSystems/RCROutflow.cpp
       ../EquationSystems/UndefinedInOutflow.h
       ../EquationSystems/UndefinedInOutflow.cpp)


ADD_SOLVER_EXECUTABLE(Fld2Tecplot solvers ${Fld2TecplotSources})

