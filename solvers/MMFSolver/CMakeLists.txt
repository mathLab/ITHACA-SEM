CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_MMF
    "Build the ADR Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_MMF )
    SET(MMFSolverSource 
       ./MMFSolver.cpp 
       ./EquationSystems/MMFMaxwell.cpp)		

    ADD_SOLVER_EXECUTABLE(MMFSolver solvers 
        ${MMFSolverSource})

ADD_NEKTAR_TEST        (TestMaxwellSphere)

ENDIF( NEKTAR_SOLVER_MMF )
