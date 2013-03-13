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
       ./RiemannSolvers/AverageSolver.cpp
       ./RiemannSolvers/AUSM0Solver.cpp
       ./RiemannSolvers/AUSM1Solver.cpp
       ./RiemannSolvers/AUSM2Solver.cpp
       ./RiemannSolvers/AUSM3Solver.cpp
       ./RiemannSolvers/CompressibleSolver.cpp
       ./RiemannSolvers/ExactSolverToro.cpp
       ./RiemannSolvers/HLLSolver.cpp
       ./RiemannSolvers/HLLCSolver.cpp
       ./RiemannSolvers/LaxFriedrichsSolver.cpp
       ./RiemannSolvers/RoeSolver.cpp
       )

    ADD_SOLVER_EXECUTABLE(CompressibleFlowSolver solvers 
			${CompressibleFlowSolverSource})

    ADD_NEKTAR_TEST(CylinderSubsonic_FRDG_SEM)
    ADD_NEKTAR_TEST(CylinderSubsonicMix)
    ADD_NEKTAR_TEST(CylinderSubsonic_P3)
    ADD_NEKTAR_TEST(CylinderSubsonic_P8)
    ADD_NEKTAR_TEST(IsentropicVortex_FRDG_SEM)
    ADD_NEKTAR_TEST(IsentropicVortex_FRSD_SEM)
    ADD_NEKTAR_TEST(IsentropicVortex_FRHU_SEM)
    ADD_NEKTAR_TEST(IsentropicVortex16_P3)
    ADD_NEKTAR_TEST(IsentropicVortex16_P8)
    ADD_NEKTAR_TEST(RinglebFlow_P3)
    ADD_NEKTAR_TEST(RinglebFlow_P8)
    ADD_NEKTAR_TEST(Couette_WeakDG_LDG_SEM)
    ADD_NEKTAR_TEST(Couette_FRHU_LDG_SEM)
    ADD_NEKTAR_TEST(Couette_WeakDG_LDG_MODIFIED)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_LDG_SEM) 
ENDIF( NEKTAR_SOLVER_COMPRESSIBLE_FLOW )
