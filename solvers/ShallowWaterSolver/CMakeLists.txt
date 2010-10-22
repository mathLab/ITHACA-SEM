SET(ShallowWaterSolverSource 
       ./ShallowWaterSolver.cpp
       ../ADRSolver/EquationSystem.cpp
       ../ADRSolver/SessionReader.cpp
       ./EquationSystems/ShallowWaterSystem.cpp
       ./EquationSystems/LinearSWE.cpp
       ./EquationSystems/NonlinearSWE.cpp		
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(ShallowWaterSolver 
			${ShallowWaterSolverSource})
