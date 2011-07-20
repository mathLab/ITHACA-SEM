SET(ShallowWaterSolverSource 
       ./ShallowWaterSolver.cpp
       ./EquationSystems/ShallowWaterSystem.cpp
       ./EquationSystems/LinearSWE.cpp
       ./EquationSystems/NonlinearSWE.cpp		
       ../Auxiliary/EquationSystem.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(ShallowWaterSolver 
			${ShallowWaterSolverSource})
