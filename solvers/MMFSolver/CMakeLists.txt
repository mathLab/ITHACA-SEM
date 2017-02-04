CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_MMF
    "Build the ADR Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_MMF )
    SET(MMFSolverSource 
       ./MMFSolver.cpp 
       ./EquationSystems/MMFAdvection.cpp 
       )		

    ADD_SOLVER_EXECUTABLE(MMFSolver solvers 
        ${MMFSolverSource})

ENDIF( NEKTAR_SOLVER_MMF )
