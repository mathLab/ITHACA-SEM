SET(ADR2DManifoldSolverSource 
       ./ADR2DManifoldSolver.cpp 
       ./ADR2DManifold.cpp
       ../Auxiliary/EquationSystem.cpp)

ADD_SOLVER_EXECUTABLE(ADR2DManifoldSolver 
			${ADR2DManifoldSolverSource})
