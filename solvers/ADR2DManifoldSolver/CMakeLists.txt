CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_ADR2DMANIFOLD
    "Build the ADR 2D manifold solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_ADR2DMANIFOLD )
    SET(ADR2DManifoldSolverSource 
       ./ADR2DManifoldSolver.cpp 
       ./ADR2DManifold.cpp)

    ADD_SOLVER_EXECUTABLE(ADR2DManifoldSolver solvers-extra 
			${ADR2DManifoldSolverSource})

ENDIF( NEKTAR_SOLVER_ADR2DMANIFOLD )
