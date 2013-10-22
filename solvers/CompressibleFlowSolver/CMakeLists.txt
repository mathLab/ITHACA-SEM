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

    ADD_NEKTAR_TEST        (CylinderSubsonic_FRDG_SEM)
    ADD_NEKTAR_TEST        (CylinderSubsonic_FRDG_GAUSS)
    ADD_NEKTAR_TEST_LENGTHY(CylinderSubsonic_FRHU_SEM)
    ADD_NEKTAR_TEST        (CylinderSubsonicMix)
    ADD_NEKTAR_TEST        (CylinderSubsonic_P3)
    ADD_NEKTAR_TEST_LENGTHY(CylinderSubsonic_P8)
    ADD_NEKTAR_TEST        (IsentropicVortex16_P3)
    ADD_NEKTAR_TEST        (IsentropicVortex_FRDG_SEM)
    ADD_NEKTAR_TEST_LENGTHY(IsentropicVortex_FRSD_SEM)
    ADD_NEKTAR_TEST_LENGTHY(IsentropicVortex_FRHU_SEM)
    ADD_NEKTAR_TEST_LENGTHY(IsentropicVortex16_P8)
    ADD_NEKTAR_TEST_LENGTHY(IsentropicVortex16_P3_GAUSS)
    ADD_NEKTAR_TEST_LENGTHY(IsentropicVortex16_P8_GAUSS)
    ADD_NEKTAR_TEST        (RinglebFlow_P3)
    ADD_NEKTAR_TEST_LENGTHY(RinglebFlow_P8)
    #ADD_NEKTAR_TEST        (Couette_WeakDG_LDG_MODIFIED)
    ADD_NEKTAR_TEST        (Couette_WeakDG_LDG_SEM)
    #ADD_NEKTAR_TEST        (Couette_WeakDG_LDG_GAUSS)
    #ADD_NEKTAR_TEST        (Couette_FRDG_LFRDG_GAUSS)
    #ADD_NEKTAR_TEST_LENGTHY(Couette_FRDG_LDG_GAUSS)
    ADD_NEKTAR_TEST_LENGTHY(Couette_FRHU_LDG_SEM)
    ADD_NEKTAR_TEST        (Couette_WeakDG_LFRHU_SEM)
    ADD_NEKTAR_TEST_LENGTHY(Couette_FRHU_LFRHU_SEM)
    ADD_NEKTAR_TEST_LENGTHY(CylinderSubsonic_NS_WeakDG_LDG_SEM)
    ADD_NEKTAR_TEST_LENGTHY(CylinderSubsonic_NS_WeakDG_LDG_GAUSS) 
    ADD_NEKTAR_TEST_LENGTHY(CylinderSubsonic_NS_FRDG_LFRDG_GAUSS)   
ENDIF( NEKTAR_SOLVER_COMPRESSIBLE_FLOW )
