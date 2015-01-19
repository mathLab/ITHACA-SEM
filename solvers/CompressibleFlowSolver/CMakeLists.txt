CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_COMPRESSIBLE_FLOW
    "Build the Compressible Flow Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_COMPRESSIBLE_FLOW )
    SET(CompressibleFlowSolverSource 
       ./CompressibleFlowSolver.cpp 
       ./EquationSystems/CompressibleFlowSystem.cpp
       ./EquationSystems/EulerCFE.cpp
       ./EquationSystems/EulerADCFE.cpp
       ./EquationSystems/NavierStokesCFE.cpp
       ./Filters/FilterEnergy.cpp
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

    SUBDIRS(Utilities)
   
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
    ADD_NEKTAR_TEST        (IsentropicVortex_WeakDG_HexDeformed)
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
    ADD_NEKTAR_TEST        (CylinderSubsonic_NS_WeakDG_LDG_variableMu)
    ADD_NEKTAR_TEST        (Couette_WeakDG_LDG_SEM_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST        (CylinderSubsonic_NS_WeakDG_LDG_SEM_3DHomo1D_MVM)
    ADD_NEKTAR_TEST        (Couette_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST        (CylinderSubsonic_NS_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST        (Couette_FRHU_LFRHU_SEM_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST        (CylinderSubsonic_NS_FRHU_LFRHU_SEM_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST        (Couette_FRSD_LFRSD_MODIFIED_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST        (CylinderSubsonic_NS_FRSD_LFRSD_MODIFIED_3DHOMO1D_MVM)

    IF (NEKTAR_USE_MPI)
        #ADD_NEKTAR_TEST_LENGTHY(Perturbation_M05_square_CBC_par)
        #ADD_NEKTAR_TEST_LENGTHY(Perturbation_M05_square_CBC_back_par)
        #ADD_NEKTAR_TEST_LENGTHY(Perturbation_M15_square_CBC_par)
        #ADD_NEKTAR_TEST_LENGTHY(Perturbation_M15_square_CBC_back_par)
        ADD_NEKTAR_TEST_LENGTHY(Perturbation_M05_circle_CBC_back_par)
        ADD_NEKTAR_TEST_LENGTHY(Perturbation_M15_circle_CBC_back_par)
        #ADD_NEKTAR_TEST_LENGTHY(Perturbation_NS_M05_square_CBC_par)
        #ADD_NEKTAR_TEST_LENGTHY(Perturbation_NS_M15_square_CBC_par)
        ADD_NEKTAR_TEST(RAE5240_BSF_LE_bcsfromfile_par)
	ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_LDG_SEM_VariableMu_par)
	ADD_NEKTAR_TEST(Couette_FRDG_LDG_GAUSS_par)
	ADD_NEKTAR_TEST(Couette_WeakDG_LDG_SEM_VariableMu_par)
	ADD_NEKTAR_TEST(Couette_FRHU_LDG_SEM_par)
    ENDIF (NEKTAR_USE_MPI)

    IF (NEKTAR_USE_FFTW)
	ADD_NEKTAR_TEST(Couette_WeakDG_LDG_SEM_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_LDG_SEM_3DHomo1D_FFT)
	ADD_NEKTAR_TEST(IsentropicVortex16_WeakDG_SEM_3DHomo1D_FFT)
        ADD_NEKTAR_TEST(Couette_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(IsentropicVortex16_FRDG_GLL_LAGRANGE_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(Couette_FRHU_LFRHU_SEM_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRHU_LFRHU_SEM_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(IsentropicVortex16_FRHU_SEM_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(Couette_FRSD_LFRSD_MODIFIED_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRSD_LFRSD_MODIFIED_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(IsentropicVortex16_FRSD_MODIFIED_3DHOMO1D_FFT)
    ENDIF (NEKTAR_USE_FFTW)
ENDIF( NEKTAR_SOLVER_COMPRESSIBLE_FLOW )
