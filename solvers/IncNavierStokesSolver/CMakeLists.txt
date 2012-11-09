SET(IncNavierStokesSolverSource 
       ./EquationSystems/CFLtester.cpp
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./AdvectionTerms/AdjointAdvection.cpp
       ./AdvectionTerms/AdvectionTerm.cpp
       ./AdvectionTerms/LinearisedAdvection.cpp
       ./AdvectionTerms/NavierStokesAdvection.cpp
       ./AdvectionTerms/SkewSymmetricAdvection.cpp
       ./AdvectionTerms/NoAdvection.cpp
       ./IncNavierStokesSolver.cpp 
       )

ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver solvers 
			${IncNavierStokesSolverSource})

ADD_NEKTAR_TEST(Channel_Flow_3modes)
ADD_NEKTAR_TEST(Channel_Flow_8modes_Body_Force)
ADD_NEKTAR_TEST(Channel_Flow_8modes)
ADD_NEKTAR_TEST(Channel_Flow_2D_BC_from_file)
ADD_NEKTAR_TEST(Kovasznay_Flow_3modes)
ADD_NEKTAR_TEST(ChanStability)
#ADD_NEKTAR_TEST(ChanStability_adj) ####
#ADD_NEKTAR_TEST(ChanStability_Coupled_3D) ###
ADD_NEKTAR_TEST(Test_2DFlow_lineforcing_bcfromfile)
ADD_NEKTAR_TEST(Test_ChanFlow2D_bcsfromfiles)
ADD_NEKTAR_TEST(Test_ChanFlow_3DH1D_MVM)
ADD_NEKTAR_TEST(Test_ChanFlow_3DH1D_Parallel_mode1)
#ADD_NEKTAR_TEST()
#ADD_NEKTAR_TEST()
#ADD_NEKTAR_TEST()
#ADD_NEKTAR_TEST()
#ADD_NEKTAR_TEST()
#ADD_NEKTAR_TEST()
#ADD_NEKTAR_TEST()



#IF (NEKTAR_USING_ARPACK)
#ADD_NEKTAR_TEST(ChanStability_adj_Ar) ###
#ADD_NEKTAR_TEST(ChanStability_Coupled) ###
#ENDIF (NEKTAR_USING_ARPACK)

IF (NEKTAR_USING_FFTW)
    ADD_NEKTAR_TEST(Test_ChanFlow_3DH1D_FFT)
    ADD_NEKTAR_TEST(Test_ChanFlow_3DH2D_FFT)
ENDIF (NEKTAR_USING_FFTW)

IF (NEKTAR_USE_MPI)
    ADD_NEKTAR_TEST(Channel_Flow_3modes_Parallel)
ENDIF (NEKTAR_USE_MPI)
