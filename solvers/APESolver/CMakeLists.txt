SET(APESolverSource 
       ./APESolver.cpp
       ./EquationSystems/APESystem.cpp
       ./EquationSystems/APE.cpp		
       ../Auxiliary/EquationSystem.cpp)

ADD_SOLVER_EXECUTABLE(APESolver solvers 
			${APESolverSource})
