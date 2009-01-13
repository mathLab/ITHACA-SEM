SET(ShallowWaterSolverSource 
       ./ShallowWaterSolver.cpp 
       ./ShallowWaterEquations.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(ShallowWaterSolver 
			${ShallowWaterSolverSource})
