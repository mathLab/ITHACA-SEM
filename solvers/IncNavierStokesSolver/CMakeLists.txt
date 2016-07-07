CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_INCNAVIERSTOKES
    "Build the Incompressible Navier-Stokes solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_INCNAVIERSTOKES )
    SET(IncNavierStokesSolverSource 
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./EquationSystems/VCSMapping.cpp
       ./EquationSystems/Extrapolate.cpp
       ./EquationSystems/StandardExtrapolate.cpp
       ./EquationSystems/MappingExtrapolate.cpp
       ./EquationSystems/SubSteppingExtrapolate.cpp
       ./AdvectionTerms/AdjointAdvection.cpp
       ./AdvectionTerms/LinearisedAdvection.cpp
       ./AdvectionTerms/NavierStokesAdvection.cpp
       ./AdvectionTerms/SkewSymmetricAdvection.cpp
       ./AdvectionTerms/AlternateSkewAdvection.cpp
       ./AdvectionTerms/NoAdvection.cpp
       ./Filters/FilterEnergy.cpp
       ./Filters/FilterReynoldsStresses.cpp
       ./Filters/FilterMovingBody.cpp
       ./Forcing/ForcingMovingBody.cpp	
       ./Forcing/ForcingStabilityCoupledLNS.cpp	
       ./IncNavierStokesSolver.cpp 
       )

    ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver solvers 
			${IncNavierStokesSolverSource})

    SUBDIRS(Utilities)      

    ADD_NEKTAR_TEST(Kovasznay_Flow_3modes)
    ADD_NEKTAR_TEST(ChanStability)
    #ADD_NEKTAR_TEST(ChanStability_adj)
    #ADD_NEKTAR_TEST(ChanStability_Coupled_3D)
    ADD_NEKTAR_TEST(2DFlow_lineforcing_bcfromfile)
    ADD_NEKTAR_TEST(ChanFlow2D_bcsfromfiles)
    ADD_NEKTAR_TEST_LENGTHY(ChanFlow_3DH1D_MVM)
    ADD_NEKTAR_TEST_LENGTHY(ChanFlow_3DH2D_MVM)
    ADD_NEKTAR_TEST(ChanFlow_LinNS_m8)
    ADD_NEKTAR_TEST(ChanFlow_m3)
    ADD_NEKTAR_TEST(ChanFlow_m3_SKS)
    ADD_NEKTAR_TEST(ChanFlow_m8)
    ADD_NEKTAR_TEST(ChanFlow_m8_BodyForce)
    ADD_NEKTAR_TEST(ChanFlow_m8_singular)
    ADD_NEKTAR_TEST(Channel_Flow_3modes_rad)    
    ADD_NEKTAR_TEST(Couette_3DH2D_MVM)
    ADD_NEKTAR_TEST(Hex_channel_m3)
    ADD_NEKTAR_TEST(Pyr_channel_m3)
    ADD_NEKTAR_TEST_LENGTHY(Hex_channel_m8)
    ADD_NEKTAR_TEST_LENGTHY(Hex_channel_m8_srhs)
    ADD_NEKTAR_TEST(Hex_Kovasnay_SubStep)
    ADD_NEKTAR_TEST(HM)
    ADD_NEKTAR_TEST(HM_Adj)
    ADD_NEKTAR_TEST_LENGTHY(bfs_tg)
    ADD_NEKTAR_TEST(SM)
    ADD_NEKTAR_TEST(SM_base_file)
    ADD_NEKTAR_TEST(SM_Adj)
    ADD_NEKTAR_TEST(SM_Adj_base_file)
    ADD_NEKTAR_TEST(KovaFlow_m3)
    ADD_NEKTAR_TEST(KovaFlow_m8)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_HOBC)
    ADD_NEKTAR_TEST(KovaFlow_varP)
    #ADD_NEKTAR_TEST(KovaFlow_Oseen_m8)
    ADD_NEKTAR_TEST_LENGTHY(KovaFlow_3DH1D_P5_20modes_MVM)
    ADD_NEKTAR_TEST_LENGTHY(KovaFlow_3DH1D_P5_20modes_MVM_Deal)
    ADD_NEKTAR_TEST_LENGTHY(KovaFlow_3DH1D_P5_20modes_SKS_MVM)
    ADD_NEKTAR_TEST(KovaFlow_SubStep_2order)
    ADD_NEKTAR_TEST(Kovas_Quad6_Tri4_mixedbcs)
    ADD_NEKTAR_TEST(SinCos_LinNS_3DHom1D)
    ADD_NEKTAR_TEST(TaylorVor_dt1)
    ADD_NEKTAR_TEST_LENGTHY(TaylorVor_dt2)
    ADD_NEKTAR_TEST(Tet_Kovasnay_SubStep)
    ADD_NEKTAR_TEST(Tet_Kovasnay_HOBC)
    ADD_NEKTAR_TEST(Tet_channel_m3)
    ADD_NEKTAR_TEST(Tet_channel_m3_dealiasing)
    ADD_NEKTAR_TEST(Tet_channel_robin_m3)
    ADD_NEKTAR_TEST(Tet_channel_m4_per)
    ADD_NEKTAR_TEST_LENGTHY(Tet_prism_channel_rad)
    ADD_NEKTAR_TEST_LENGTHY(Tet_equitri)
    ADD_NEKTAR_TEST_LENGTHY(Prism_channel_m6)
    ADD_NEKTAR_TEST_LENGTHY(Prism_channel_m6_dealiasing)
    ADD_NEKTAR_TEST(ChannelSpongeNSE)
    ADD_NEKTAR_TEST(ChannelSpongeLNSE)
    ADD_NEKTAR_TEST(ChanFlow_Standard_BodyForce)
    ADD_NEKTAR_TEST(Cyl_AdaptiveSFD)
    ADD_NEKTAR_TEST(CylFlow_Mov_mapping)
    ADD_NEKTAR_TEST(PPF_R10000_ModifiedArnoldi_Shift) 
    ADD_NEKTAR_TEST(PPF_R15000_ModifiedArnoldi_Shift) 

    IF (NEKTAR_USE_ARPACK)
        ADD_NEKTAR_TEST(ChanStability_adj_Ar) 
        ADD_NEKTAR_TEST(bfs_tg-AR)
        ADD_NEKTAR_TEST(ChanStability_Coupled) 
        ADD_NEKTAR_TEST(PPF_R10000_Arpack_LM) 
        ADD_NEKTAR_TEST(PPF_R15000_Arpack_Shift) 
        #ADD_NEKTAR_TEST(PPF_R15000_Arpack_NoImagShift) 
        ADD_NEKTAR_TEST(PPF_R15000_Arpack_NoImagShift_LM) 
    ENDIF (NEKTAR_USE_ARPACK)

    IF (NEKTAR_USE_FFTW)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_FFT)
        ADD_NEKTAR_TEST(ChanFlow_3DH2D_FFT)
        ADD_NEKTAR_TEST(CylFlow_MovBody)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_MVM_FFTW_Consistency)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_P8_16modes_Mapping-explicit)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_P8_16modes_Mapping-implicit)
        ADD_NEKTAR_TEST(KovaFlow_m8_short_HOBC_3D1H)
        ADD_NEKTAR_TEST_LENGTHY(KovaFlow_3DH1D_adaptive_16modes_FFTW_Mapping)
        ADD_NEKTAR_TEST(Couette_3DH2D_FFT)
    ENDIF (NEKTAR_USE_FFTW)

    IF (NEKTAR_USE_MPI)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Parallel_mode1)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Parallel_mode2)
        ADD_NEKTAR_TEST(ChanFlow_m3_par)
        ADD_NEKTAR_TEST_LENGTHY(ChanFlow_m8_BodyForce_par)
        ADD_NEKTAR_TEST_LENGTHY(Hex_channel_m8_par)
        ADD_NEKTAR_TEST_LENGTHY(Pyr_channel_m6_par)
        ADD_NEKTAR_TEST_LENGTHY(Hex_channel_m8_srhs_par)
        ADD_NEKTAR_TEST(Tet_channel_m4_per_xxt_ml_par)
        ADD_NEKTAR_TEST(Tet_channel_m8_par)
        ADD_NEKTAR_TEST_LENGTHY(Tet_channel_m8_iter_ml_par)
        ADD_NEKTAR_TEST_LENGTHY(bfs_tg_par)
        IF (NEKTAR_USE_PETSC)
            ADD_NEKTAR_TEST_LENGTHY(Tet_channel_m8_petsc_sc_par)
        ENDIF (NEKTAR_USE_PETSC)
    ENDIF (NEKTAR_USE_MPI)
ENDIF( NEKTAR_SOLVER_INCNAVIERSTOKES )
