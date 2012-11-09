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

ADD_NEKTAR_TEST(Test_Advection1D_FRDG_GLL_LAGRANGE_SEM)
ADD_NEKTAR_TEST(Test_Advection1D_FRDG_GLL_LAGRANGE)
ADD_NEKTAR_TEST(Test_Advection1D_FRDG_MODIFIED)
ADD_NEKTAR_TEST(Test_Advection1D_FRSD_GLL_LAGRANGE_SEM)
ADD_NEKTAR_TEST(Test_Advection1D_FRSD_GLL_LAGRANGE)
ADD_NEKTAR_TEST(Test_Advection1D_FRSD_MODIFIED)
ADD_NEKTAR_TEST(Test_Advection1D_FRHU_GLL_LAGRANGE_SEM)
ADD_NEKTAR_TEST(Test_Advection1D_FRHU_GLL_LAGRANGE)
ADD_NEKTAR_TEST(Test_Advection1D_FRHU_MODIFIED)
ADD_NEKTAR_TEST(Test_Advection1D_WeakDG_GLL_LAGRANGE)
ADD_NEKTAR_TEST(Test_Advection1D_WeakDG_MODIFIED)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_FRDG_GLL_LAGRANGE_SEM)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_FRDG_GLL_LAGRANGE)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_FRDG_MODIFIED)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_FRSD_GLL_LAGRANGE_SEM)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_FRSD_GLL_LAGRANGE)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_FRSD_MODIFIED)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_FRHU_GLL_LAGRANGE_SEM)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_FRHU_GLL_LAGRANGE)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_FRHU_MODIFIED)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_WeakDG_GLL_LAGRANGE)
ADD_NEKTAR_TEST(Test_InviscidBurger1D_WeakDG_MODIFIED)


ADD_NEKTAR_TEST(Test_Helmholtz2D_modal_DG)
ADD_NEKTAR_TEST(Test_Helmholtz2D_modal)
ADD_NEKTAR_TEST(Test_Helmholtz2D_nodal)
ADD_NEKTAR_TEST(Test_Helmholtz_3DHomo1D_FFT)
ADD_NEKTAR_TEST(Test_Helmholtz_3DHomo1D_MVM)
ADD_NEKTAR_TEST(Test_Helmholtz_3DHomo2D_FFT)
ADD_NEKTAR_TEST(Test_Helmholtz_3DHomo2D_MVM)
ADD_NEKTAR_TEST(Test_Helmholtz3D_modal)
ADD_NEKTAR_TEST(Test_Helmholtz3D_nodal)
ADD_NEKTAR_TEST(Test_ImDiffusion_m12)
ADD_NEKTAR_TEST(Test_ImDiffusion_m6)
ADD_NEKTAR_TEST(Test_ImDiffusion_VarCoeff)


