CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_DIFFUSION
    "Build the Diffusion Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_DIFFUSION )
    SET(DiffusionSolverSource 
        ./DiffusionSolver.cpp)
    SET(DiffusionSolverTimeIntSource
        ./DiffusionSolverTimeInt.cpp)
        
    ADD_SOLVER_EXECUTABLE(DiffusionSolver solvers 
        ${DiffusionSolverSource})
    ADD_SOLVER_EXECUTABLE(DiffusionSolverTimeInt solvers 
        ${DiffusionSolverTimeIntSource})

    ADD_NEKTAR_TEST(ImDiffusion_m6)
    ADD_NEKTAR_TEST(ImDiffusion_m6_time_int)
ENDIF( NEKTAR_SOLVER_DIFFUSION )
