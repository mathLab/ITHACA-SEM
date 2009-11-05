SET(ADR2DManifoldSolverSource 
       ./ADR2DManifoldSolver.cpp 
       ./ADR2DManifold.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(ADR2DManifoldSolver 
			${ADR2DManifoldSolverSource})
