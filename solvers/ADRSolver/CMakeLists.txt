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

