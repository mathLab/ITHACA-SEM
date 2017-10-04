CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_DIFFUSION
    "Build the Diffusion Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_DIFFUSION OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_DIFFUSION )
    SET(DiffusionSolverSource 
        ./DiffusionSolver.cpp)
    SET(DiffusionSolverTimeIntSource
        ./DiffusionSolverTimeInt.cpp)
        
    ADD_SOLVER_EXECUTABLE(DiffusionSolver SOURCES DiffusionSolver.cpp)
    ADD_SOLVER_EXECUTABLE(DiffusionSolverTimeInt SOURCES DiffusionSolverTimeInt COMPONENT diffusionsolver)

    ADD_NEKTAR_TEST(ImDiffusion_m6)
    ADD_NEKTAR_TEST(ImDiffusion_m6_time_int)
ENDIF( NEKTAR_SOLVER_DIFFUSION )
