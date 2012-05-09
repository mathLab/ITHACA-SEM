SET(PulseWaveSolverSource 
       ./PulseWaveSolver.cpp 
       ./EquationSystems/PulseWaveSystem.cpp
       ./EquationSystems/PulseWavePropagation.cpp)

ADD_SOLVER_EXECUTABLE(PulseWaveSolver solvers 
			${PulseWaveSolverSource})
