SET(ShallowWaterSolverSource 
       ./ShallowWaterSolver.cpp
       ../ADRSolver/EquationSystem.cpp
       ./EquationSystems/ShallowWaterSystem.cpp
       ./EquationSystems/LinearSWE.cpp
       ./EquationSystems/NonlinearSWE.cpp		
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(ShallowWaterSolver 
			${ShallowWaterSolverSource})
