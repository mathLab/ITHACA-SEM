SET(ADR2DManifoldSolverSource 
       ./ADR2DManifoldSolver.cpp 
       ./ADR2DManifold.cpp)

ADD_SOLVER_EXECUTABLE(ADR2DManifoldSolver solvers-extra 
			${ADR2DManifoldSolverSource})
