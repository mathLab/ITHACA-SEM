SET(PulseWaveSolverSource 
       ./PulseWaveSolver.cpp 
       ./EquationSystems/PulseWaveSystem.cpp
       ./EquationSystems/PulseWavePropagation.cpp
       ../Auxiliary/Driver.cpp
       ../Auxiliary/DriverStandard.cpp
       ../Auxiliary/DriverArnoldi.cpp
       ../Auxiliary/DriverModifiedArnoldi.cpp
       ../Auxiliary/Filters/Filter.cpp
       ../Auxiliary/Filters/FilterThresholdMax.cpp
       ../Auxiliary/UnsteadySystem.cpp
       ../Auxiliary/EquationSystem.cpp)

IF (NEKTAR_USE_ARPACK)
    SET(PulseWaveSolverSource ${PulseWaveSolverSource}
               ../Auxiliary/DriverArpack.cpp)
ENDIF (NEKTAR_USE_ARPACK)

ADD_SOLVER_EXECUTABLE(PulseWaveSolver solvers 
			${PulseWaveSolverSource})
