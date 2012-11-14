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
#ADD_NEKTAR_TEST(ChanStability_adj)
#ADD_NEKTAR_TEST(ChanStability_Coupled_3D)
ADD_NEKTAR_TEST(Test_2DFlow_lineforcing_bcfromfile)
ADD_NEKTAR_TEST(Test_ChanFlow2D_bcsfromfiles)
ADD_NEKTAR_TEST(Test_ChanFlow_3DH1D_MVM)
ADD_NEKTAR_TEST(Test_ChanFlow_3DH1D_Parallel_mode1)
ADD_NEKTAR_TEST(Test_ChanFlow_3DH1D_Parallel_mode2)
ADD_NEKTAR_TEST(Test_ChanFlow_3DH2D_MVM)
ADD_NEKTAR_TEST(Test_ChanFlow_LinNS_m8)
ADD_NEKTAR_TEST(Test_ChanFlow_m3)
ADD_NEKTAR_TEST(Test_ChanFlow_m3_SKS)
ADD_NEKTAR_TEST(Test_ChanFlow_m8)
ADD_NEKTAR_TEST(Test_ChanFlow_m8_BodyForce)
ADD_NEKTAR_TEST(Test_ChanFlow_m8_singular)
ADD_NEKTAR_TEST(Test_Hex_Kovasnay_SubStep)
ADD_NEKTAR_TEST(Test_HM)
ADD_NEKTAR_TEST(Test_HM_Adj)
ADD_NEKTAR_TEST(bfs_tg)
ADD_NEKTAR_TEST(Test_SM)
ADD_NEKTAR_TEST(Test_SM_Adj)
ADD_NEKTAR_TEST(Test_KovaFlow_m3)
ADD_NEKTAR_TEST(Test_KovaFlow_m8)
#ADD_NEKTAR_TEST(Test_KovaFlow_Oseen_m8)
ADD_NEKTAR_TEST(Test_KovaFlow_3DH1D_P5_20modes_MVM)
ADD_NEKTAR_TEST(Test_KovaFlow_3DH1D_P5_20modes_SKS_MVM)
ADD_NEKTAR_TEST(Test_KovaFlow_SubStep_2order)
ADD_NEKTAR_TEST(Test_Kovas_Quad6_Tri4_mixedbcs)
ADD_NEKTAR_TEST(Test_SinCos_LinNS_3DHom1D)
ADD_NEKTAR_TEST(Test_TaylorVor_dt1)
ADD_NEKTAR_TEST(Test_TaylorVor_dt2)
ADD_NEKTAR_TEST(Test_Hex_channel_m3)
ADD_NEKTAR_TEST(Test_Hex_channel_m8)
ADD_NEKTAR_TEST(Test_Tet_channel_m3)
ADD_NEKTAR_TEST(Test_Tet_equitri)
ADD_NEKTAR_TEST(Test_Prism_channel_m6)



#IF (NEKTAR_USING_ARPACK)
    #ADD_NEKTAR_TEST(ChanStability_adj_Ar) ###
    #ADD_NEKTAR_TEST(ChanStability_Coupled) ###
    #ADD_NEKTAR_TEST(bfs_tg-AR)
#ENDIF (NEKTAR_USING_ARPACK)

IF (NEKTAR_USING_FFTW)
    ADD_NEKTAR_TEST(Test_ChanFlow_3DH1D_FFT)
    ADD_NEKTAR_TEST(Test_ChanFlow_3DH2D_FFT)
ENDIF (NEKTAR_USING_FFTW)

IF (NEKTAR_USE_MPI)
    ADD_NEKTAR_TEST(Channel_Flow_3modes_Parallel)
    ADD_NEKTAR_TEST(Test_ChanFlow_m3_par)
    ADD_NEKTAR_TEST(Test_Hex_channel_m8_par)
    ADD_NEKTAR_TEST(Test_Tet_channel_m8_par)
ENDIF (NEKTAR_USE_MPI)
