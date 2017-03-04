CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_DIFFUSION
    "Build the Diffusion Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_DIFFUSION )
    SET(DiffusionSolverSource 
        ./DiffusionSolver.cpp)
    SET(DiffusionSolverTimeIntSource
        ./DiffusionSolverTimeInt.cpp)
    SET(MMFDiffusionSource 
        ./MMFDiffusion.cpp)
        
    ADD_SOLVER_EXECUTABLE(DiffusionSolver solvers 
        ${DiffusionSolverSource})
    ADD_SOLVER_EXECUTABLE(DiffusionSolverTimeInt solvers 
        ${DiffusionSolverTimeIntSource})
    ADD_SOLVER_EXECUTABLE(MMFDiffusion solvers 
        ${MMFDiffusionSource})

    #MMF Diffusion
    ADD_NEKTAR_TEST(TestMMFDiffPlane)
    ADD_NEKTAR_TEST(TestMMFDiffSphere)
    ADD_NEKTAR_TEST(TestMMFDiffCube)

    ADD_NEKTAR_TEST(ImDiffusion_m6)
    ADD_NEKTAR_TEST(ImDiffusion_m6_time_int)
ENDIF( NEKTAR_SOLVER_DIFFUSION )
