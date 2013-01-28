CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_COMPRESSIBLE_FLOW
    "Build the Compressible Flow Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_COMPRESSIBLE_FLOW )
    SET(CompressibleFlowSolverSource 
       ./CompressibleFlowSolver.cpp 
       ./EquationSystems/CompressibleFlowSystem.cpp
       ./EquationSystems/EulerCFE.cpp
       ./EquationSystems/EulerArtificialDiffusionCFE.cpp
       ./EquationSystems/NavierStokesCFE.cpp
       ./RiemannSolvers/CompressibleSolver.cpp
       ./RiemannSolvers/AverageSolver.cpp
       ./RiemannSolvers/ExactSolver.cpp)

    ADD_SOLVER_EXECUTABLE(CompressibleFlowSolver solvers 
			${CompressibleFlowSolverSource})

    ADD_NEKTAR_TEST(CylinderSubsonicMix)
    ADD_NEKTAR_TEST(CylinderSubsonic_P3)
    ADD_NEKTAR_TEST(CylinderSubsonic_P8)
    ADD_NEKTAR_TEST(IsentropicVortex16_P3)
    ADD_NEKTAR_TEST(IsentropicVortex16_P8)
    ADD_NEKTAR_TEST(RinglebFlow_P3)
    ADD_NEKTAR_TEST(RinglebFlow_P8)
ENDIF( NEKTAR_SOLVER_COMPRESSIBLE_FLOW )
