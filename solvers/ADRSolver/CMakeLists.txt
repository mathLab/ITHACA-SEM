CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_ADR
    "Build the ADR Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_ADR )
    SET(ADRSolverSource 
       ./ADRSolver.cpp 
       ./EquationSystems/Laplace.cpp
       ./EquationSystems/Poisson.cpp       
       ./EquationSystems/Helmholtz.cpp
       ./EquationSystems/SteadyAdvectionDiffusion.cpp
       ./EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ./EquationSystems/UnsteadyAdvection.cpp
       ./EquationSystems/UnsteadyDiffusion.cpp
       ./EquationSystems/UnsteadyAdvectionDiffusion.cpp
       ./EquationSystems/UnsteadyInviscidBurger.cpp
       ./EquationSystems/CFLtester.cpp
       ./EquationSystems/EigenValuesAdvection.cpp)

    ADD_SOLVER_EXECUTABLE(ADRSolver solvers 
			${ADRSolverSource})

    ADD_NEKTAR_TEST(Advection1D_FRDG_GLL_LAGRANGE_SEM)
    ADD_NEKTAR_TEST(Advection1D_FRDG_GLL_LAGRANGE)
    ADD_NEKTAR_TEST(Advection1D_FRDG_MODIFIED)
    ADD_NEKTAR_TEST(Advection1D_FRSD_GLL_LAGRANGE_SEM)
    ADD_NEKTAR_TEST(Advection1D_FRSD_GLL_LAGRANGE)
    ADD_NEKTAR_TEST(Advection1D_FRSD_MODIFIED)
    ADD_NEKTAR_TEST(Advection1D_FRHU_GLL_LAGRANGE_SEM)
    ADD_NEKTAR_TEST(Advection1D_FRHU_GLL_LAGRANGE)
    ADD_NEKTAR_TEST(Advection1D_FRHU_MODIFIED)
    ADD_NEKTAR_TEST(Advection1D_WeakDG_GLL_LAGRANGE)
    ADD_NEKTAR_TEST(Advection1D_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(InviscidBurger1D_FRDG_GLL_LAGRANGE_SEM)
    ADD_NEKTAR_TEST(InviscidBurger1D_FRDG_GLL_LAGRANGE)
    ADD_NEKTAR_TEST(InviscidBurger1D_FRDG_MODIFIED)
    ADD_NEKTAR_TEST(InviscidBurger1D_FRSD_GLL_LAGRANGE_SEM)
    ADD_NEKTAR_TEST(InviscidBurger1D_FRSD_GLL_LAGRANGE)
    ADD_NEKTAR_TEST(InviscidBurger1D_FRSD_MODIFIED)
    ADD_NEKTAR_TEST(InviscidBurger1D_FRHU_GLL_LAGRANGE_SEM)
    ADD_NEKTAR_TEST(InviscidBurger1D_FRHU_GLL_LAGRANGE)
    ADD_NEKTAR_TEST(InviscidBurger1D_FRHU_MODIFIED)
    ADD_NEKTAR_TEST(InviscidBurger1D_WeakDG_GLL_LAGRANGE)
    ADD_NEKTAR_TEST(InviscidBurger1D_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(Advection2D_dirichlet_deformed_GLL_LAGRANGE_10x10)
    ADD_NEKTAR_TEST(Advection2D_dirichlet_deformed_GLL_LAGRANGE_SEM_10x10)
    ADD_NEKTAR_TEST(Advection2D_dirichlet_deformed_MODIFIED_10x10)
    ADD_NEKTAR_TEST(Advection2D_dirichlet_regular_GLL_LAGRANGE_10x10)
    ADD_NEKTAR_TEST(Advection2D_dirichlet_regular_GLL_LAGRANGE_SEM_10x10)
    ADD_NEKTAR_TEST(Advection2D_dirichlet_regular_MODIFIED_10x10)
    ADD_NEKTAR_TEST(Advection2D_ISO_deformed_GLL_LAGRANGE_3x3)
    ADD_NEKTAR_TEST(Advection2D_ISO_deformed_GLL_LAGRANGE_SEM_3x3)
    ADD_NEKTAR_TEST(Advection2D_ISO_deformed_MODIFIED_3x3)
    ADD_NEKTAR_TEST(Advection2D_ISO_regular_GLL_LAGRANGE_3x3)
    ADD_NEKTAR_TEST(Advection2D_ISO_regular_GLL_LAGRANGE_SEM_3x3)
    ADD_NEKTAR_TEST(Advection2D_ISO_regular_MODIFIED_3x3)
    ADD_NEKTAR_TEST(Advection2D_periodic_deformed_GLL_LAGRANGE_2x2)
    ADD_NEKTAR_TEST(Advection2D_periodic_deformed_GLL_LAGRANGE_SEM_2x2)
    ADD_NEKTAR_TEST(Advection2D_periodic_deformed_MODIFIED_2x2)
    ADD_NEKTAR_TEST(Advection2D_periodic_regular_GLL_LAGRANGE_2x2)
    ADD_NEKTAR_TEST(Advection2D_periodic_regular_GLL_LAGRANGE_SEM_2x2)
    ADD_NEKTAR_TEST(Advection2D_periodic_regular_MODIFIED_2x2)
    ADD_NEKTAR_TEST(Advection2D_periodic_deformed_GLL_LAGRANGE_10x10)
    ADD_NEKTAR_TEST(Advection2D_periodic_deformed_GLL_LAGRANGE_SEM_10x10)
    ADD_NEKTAR_TEST(Advection2D_periodic_deformed_MODIFIED_10x10)
    ADD_NEKTAR_TEST(Advection3D_m12_DG_hex_periodic)
    ADD_NEKTAR_TEST(Advection3D_m12_DG_hex)
    ADD_NEKTAR_TEST(Advection3D_m12_DG_prism)
    ADD_NEKTAR_TEST(Advection3D_m12_DG_tet)
    ADD_NEKTAR_TEST(Advection_m12_DG_Order1)
    ADD_NEKTAR_TEST(Advection_m12_DG_Order2)
    ADD_NEKTAR_TEST(Advection_m12_DG_periodic)
    ADD_NEKTAR_TEST(Advection_m12_Order1)
    ADD_NEKTAR_TEST(Advection_m12_Order2)
    ADD_NEKTAR_TEST(Advection_m14_DG_Order4)
    ADD_NEKTAR_TEST(Advection_m14_Order4)
    ADD_NEKTAR_TEST(ExDiffusion_m3)
    ADD_NEKTAR_TEST(ExDiffusion_m8)
    ADD_NEKTAR_TEST(Helmholtz1D_8modes_DG)
    ADD_NEKTAR_TEST(Helmholtz1D_8modes)
    ADD_NEKTAR_TEST(Helmholtz1D_8nodes)
    
    ADD_NEKTAR_TEST(SteadyAdvDiffReact2D_modal)
    ADD_NEKTAR_TEST(UnsteadyAdvectionDiffusion_3DHomo1D_MVM)
    ADD_NEKTAR_TEST(UnsteadyAdvectionDiffusion_3DHomo2D_MVM)
    ADD_NEKTAR_TEST(UnsteadyAdvectionDiffusion_Order1_0001)
    ADD_NEKTAR_TEST(UnsteadyAdvectionDiffusion_Order1_001)
    ADD_NEKTAR_TEST(UnsteadyAdvectionDiffusion_Order2_0001)
    ADD_NEKTAR_TEST(UnsteadyAdvectionDiffusion_Order2_001)
    
    ADD_NEKTAR_TEST(Helmholtz2D_modal_DG)
    ADD_NEKTAR_TEST(Helmholtz2D_modal)
    ADD_NEKTAR_TEST(Helmholtz2D_nodal)
    ADD_NEKTAR_TEST(Helmholtz_3DHomo1D_MVM)
    ADD_NEKTAR_TEST(Helmholtz_3DHomo2D_MVM)
    ADD_NEKTAR_TEST(Helmholtz3D_modal)
    ADD_NEKTAR_TEST(Helmholtz3D_nodal)
    ADD_NEKTAR_TEST(ImDiffusion_m12)
    ADD_NEKTAR_TEST(ImDiffusion_m6)
    ADD_NEKTAR_TEST(ImDiffusion_VarCoeff)
    
    IF (NEKTAR_USE_FFTW)
        ADD_NEKTAR_TEST(Helmholtz_3DHomo1D_FFT)
        ADD_NEKTAR_TEST(Helmholtz_3DHomo2D_FFT)
        ADD_NEKTAR_TEST(UnsteadyAdvectionDiffusion_3DHomo1D_FFT)
        ADD_NEKTAR_TEST(UnsteadyAdvectionDiffusion_3DHomo2D_FFT)
    ENDIF (NEKTAR_USE_FFTW)

ENDIF( NEKTAR_SOLVER_ADR )
