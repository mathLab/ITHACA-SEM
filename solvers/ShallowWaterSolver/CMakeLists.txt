SET(ShallowWaterSolverSource 
       ./ShallowWaterSolver.cpp
       ./EquationSystems/ShallowWaterSystem.cpp
       ./EquationSystems/LinearSWE.cpp
       ./EquationSystems/NonlinearSWE.cpp		
       ../Auxiliary/EquationSystem.cpp)

ADD_SOLVER_EXECUTABLE(ShallowWaterSolver solvers 
			${ShallowWaterSolverSource})
