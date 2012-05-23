SET(PulseWaveSolverSource 
       ./PulseWaveSolver.cpp 
       ./EquationSystems/PulseWaveSystem.h
       ./EquationSystems/PulseWaveSystem.cpp
       ./EquationSystems/PulseWavePropagation.h
       ./EquationSystems/PulseWavePropagation.cpp)

ADD_SOLVER_EXECUTABLE(PulseWaveSolver solvers 
			${PulseWaveSolverSource})
