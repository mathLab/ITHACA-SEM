CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_INCNAVIERSTOKES
    "Build the Incompressible Navier-Stokes solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_INCNAVIERSTOKES )
    SET(IncNavierStokesSolverSource 
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
    ADD_NEKTAR_TEST(2DFlow_lineforcing_bcfromfile)
    ADD_NEKTAR_TEST(ChanFlow2D_bcsfromfiles)
    ADD_NEKTAR_TEST(ChanFlow_3DH1D_MVM)
    ADD_NEKTAR_TEST(ChanFlow_3DH2D_MVM)
    ADD_NEKTAR_TEST(ChanFlow_LinNS_m8)
    ADD_NEKTAR_TEST(ChanFlow_m3)
    ADD_NEKTAR_TEST(ChanFlow_m3_SKS)
    ADD_NEKTAR_TEST(ChanFlow_m8)
    ADD_NEKTAR_TEST(ChanFlow_m8_BodyForce)
    ADD_NEKTAR_TEST(ChanFlow_m8_singular)
    ADD_NEKTAR_TEST(Hex_Kovasnay_SubStep)
    ADD_NEKTAR_TEST(HM)
    ADD_NEKTAR_TEST(HM_Adj)
    ADD_NEKTAR_TEST(bfs_tg)
    ADD_NEKTAR_TEST(SM)
    ADD_NEKTAR_TEST(SM_Adj)
    ADD_NEKTAR_TEST(KovaFlow_m3)
    ADD_NEKTAR_TEST(KovaFlow_m8)
    #ADD_NEKTAR_TEST(KovaFlow_Oseen_m8)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_MVM)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_SKS_MVM)
    ADD_NEKTAR_TEST(KovaFlow_SubStep_2order)
    ADD_NEKTAR_TEST(Kovas_Quad6_Tri4_mixedbcs)
    ADD_NEKTAR_TEST(SinCos_LinNS_3DHom1D)
    ADD_NEKTAR_TEST(TaylorVor_dt1)
    ADD_NEKTAR_TEST(TaylorVor_dt2)
    ADD_NEKTAR_TEST(Hex_channel_m3)
    ADD_NEKTAR_TEST(Hex_channel_m8)
    ADD_NEKTAR_TEST(Hex_channel_m8_srhs)
    ADD_NEKTAR_TEST(Tet_channel_m3)
    ADD_NEKTAR_TEST(Tet_equitri)
    ADD_NEKTAR_TEST(Prism_channel_m6)

    #IF (NEKTAR_USING_ARPACK)
        #ADD_NEKTAR_TEST(ChanStability_adj_Ar) ###
        #ADD_NEKTAR_TEST(ChanStability_Coupled) ###
        #ADD_NEKTAR_TEST(bfs_tg-AR)
    #ENDIF (NEKTAR_USING_ARPACK)

    IF (NEKTAR_USING_FFTW)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_FFT)
        ADD_NEKTAR_TEST(ChanFlow_3DH2D_FFT)
    ENDIF (NEKTAR_USING_FFTW)

    IF (NEKTAR_USE_MPI)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Parallel_mode1)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Parallel_mode2)
        ADD_NEKTAR_TEST(Channel_Flow_3modes_Parallel)
        ADD_NEKTAR_TEST(ChanFlow_m3_par)
        ADD_NEKTAR_TEST(Hex_channel_m8_par)
        ADD_NEKTAR_TEST(Tet_channel_m8_par)
        ADD_NEKTAR_TEST(Tet_channel_m8_iter_ml_par)
    ENDIF (NEKTAR_USE_MPI)
ENDIF( NEKTAR_SOLVER_INCNAVIERSTOKES )
