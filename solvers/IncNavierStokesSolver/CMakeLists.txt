CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_INCNAVIERSTOKES
    "Build the Incompressible Navier-Stokes solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_INCNAVIERSTOKES OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_INCNAVIERSTOKES )
    ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver
        SOURCES ./IncNavierStokesSolver.cpp
        LIBRARY_SOURCES
        ./EquationSystems/CoupledLinearNS.cpp
        ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
        ./EquationSystems/IncNavierStokes.cpp
        ./EquationSystems/VelocityCorrectionScheme.cpp
        ./EquationSystems/VelocityCorrectionSchemeWeakPressure.cpp
        ./EquationSystems/VCSMapping.cpp
        ./EquationSystems/Extrapolate.cpp
        ./EquationSystems/StandardExtrapolate.cpp
        ./EquationSystems/MappingExtrapolate.cpp
        ./EquationSystems/SubSteppingExtrapolate.cpp
        ./EquationSystems/SubSteppingExtrapolateWeakPressure.cpp
        ./EquationSystems/WeakPressureExtrapolate.cpp
        ./AdvectionTerms/AdjointAdvection.cpp
        ./AdvectionTerms/LinearisedAdvection.cpp
        ./AdvectionTerms/NavierStokesAdvection.cpp
        ./AdvectionTerms/SkewSymmetricAdvection.cpp
        ./AdvectionTerms/AlternateSkewAdvection.cpp
        ./AdvectionTerms/NoAdvection.cpp
        ./Filters/FilterReynoldsStresses.cpp
        ./Filters/FilterMovingBody.cpp
        ./Forcing/ForcingMovingBody.cpp
        ./Forcing/ForcingStabilityCoupledLNS.cpp
        )

    SUBDIRS(Utilities)      

    ADD_NEKTAR_TEST(Kovasznay_Flow_3modes)
    ADD_NEKTAR_TEST(ChanStability)
    #ADD_NEKTAR_TEST(ChanStability_adj)
    #ADD_NEKTAR_TEST(ChanStability_Coupled_3D)
    ADD_NEKTAR_TEST(2DFlow_lineforcing_bcfromfile)
    ADD_NEKTAR_TEST(BercovierEngelman)
    ADD_NEKTAR_TEST(ChanFlow2D_bcsfromfiles)
    ADD_NEKTAR_TEST(ChanFlow_3DH1D_bcsfromfiles)
    ADD_NEKTAR_TEST(ChanFlow_3DH1D_MVM LENGTHY)
    ADD_NEKTAR_TEST(ChanFlow_3DH2D_MVM LENGTHY)
    ADD_NEKTAR_TEST(ChanFlow_LinNS_m8)
    ADD_NEKTAR_TEST(ChanFlow_m3)
    ADD_NEKTAR_TEST(ChanFlow_m3_ConOBC)
    ADD_NEKTAR_TEST(ChanFlow_m3_VCSWeakPress)
    ADD_NEKTAR_TEST(ChanFlow_m3_VCSWeakPress_ConOBC)
    ADD_NEKTAR_TEST(ChanFlow_m3_SKS)
    ADD_NEKTAR_TEST(ChanFlow_m8)
    ADD_NEKTAR_TEST(ChanFlow_m8_BodyForce)
    ADD_NEKTAR_TEST(ChanFlow_m8_singular)
    ADD_NEKTAR_TEST(ChanFlow_V8P7_Avg)
    ADD_NEKTAR_TEST(Channel_Flow_3modes_rad)    
    ADD_NEKTAR_TEST(channelTemp)    
    ADD_NEKTAR_TEST(Couette_3DH2D_MVM)
    ADD_NEKTAR_TEST(Hex_channel_m3)
    ADD_NEKTAR_TEST(Hex_channel_varP)
    ADD_NEKTAR_TEST(Pyr_channel_m3)
    ADD_NEKTAR_TEST(Pyr_channel_varP)
    ADD_NEKTAR_TEST(Hex_channel_m6_nodalRestart)
    ADD_NEKTAR_TEST(Hex_channel_m3_srhs LENGTHY)
    ADD_NEKTAR_TEST(Hex_Kovasnay_SubStep)
    ADD_NEKTAR_TEST(HM)
    ADD_NEKTAR_TEST(HM_Adj)
    ADD_NEKTAR_TEST(bfs_tg LENGTHY)
    ADD_NEKTAR_TEST(SM)
    ADD_NEKTAR_TEST(SM_base_file)
    ADD_NEKTAR_TEST(SM_Adj)
    ADD_NEKTAR_TEST(SM_Adj_base_file)
    ADD_NEKTAR_TEST(KovaFlow_m3)
    ADD_NEKTAR_TEST(KovaFlow_m8)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_ConOBC)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_ConOBC_VCSWeakPress)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_HOBC)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_HOBC_VCSWeakPress)
    ADD_NEKTAR_TEST(KovaFlow_varP)
    ADD_NEKTAR_TEST(KovaFlow_varP_per)
    #ADD_NEKTAR_TEST(KovaFlow_Oseen_m8)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_MVM LENGTHY)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_MVM_SVVHomo1D LENGTHY)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_MVM_Deal LENGTHY)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_SKS_MVM LENGTHY)
    ADD_NEKTAR_TEST(KovaFlow_SubStep_2order)
    ADD_NEKTAR_TEST(KovaFlow_SubStep_2order_VCSWeakPress)
    ADD_NEKTAR_TEST(Kovas_Quad6_Tri4_mixedbcs)
    ADD_NEKTAR_TEST(Moving_reference_frame)
    ADD_NEKTAR_TEST(SinCos_LinNS_3DHom1D)
    ADD_NEKTAR_TEST(TaylorVor_dt1)
    ADD_NEKTAR_TEST(TaylorVor_dt2 LENGTHY)
    ADD_NEKTAR_TEST(Tet_Kovasnay_SubStep)
    ADD_NEKTAR_TEST(Tet_Kovasnay_HOBC)
    ADD_NEKTAR_TEST(Tet_Kovasnay_SVV_DGKer)
    ADD_NEKTAR_TEST(Tet_channel_m3)
    ADD_NEKTAR_TEST(Tet_channel_m3_ConOBC)
    ADD_NEKTAR_TEST(Tet_channel_m3_VCSWeakPress)
    ADD_NEKTAR_TEST(Tet_channel_m3_VCSWeakPress_ConOBC)
    ADD_NEKTAR_TEST(Tet_channel_m3_dealiasing)
    ADD_NEKTAR_TEST(Tet_channel_robin_m3)
    ADD_NEKTAR_TEST(Tet_channel_m4_per)
    ADD_NEKTAR_TEST(Tet_channel_varP)
    ADD_NEKTAR_TEST(TriQuadChannel)
    ADD_NEKTAR_TEST(Tet_prism_channel_rad LENGTHY)
    ADD_NEKTAR_TEST(Tet_equitri LENGTHY)
    ADD_NEKTAR_TEST(Prism_channel_m6 LENGTHY)
    ADD_NEKTAR_TEST(Prism_channel_m6_dealiasing LENGTHY)
    ADD_NEKTAR_TEST(Prism_channel_varP)
    ADD_NEKTAR_TEST(ChannelSpongeNSE)
    ADD_NEKTAR_TEST(ChannelSpongeLNSE)
    ADD_NEKTAR_TEST(ChanFlow_Standard_BodyForce)
    ADD_NEKTAR_TEST(Cyl_AdaptiveSFD)
    ADD_NEKTAR_TEST(Womersley_PipeFlow)
    ADD_NEKTAR_TEST(CylFlow_adaptiveP LENGTHY)
    ADD_NEKTAR_TEST(CylFlow_Mov_mapping)
    ADD_NEKTAR_TEST(PPF_R10000_ModifiedArnoldi_Shift) 
    ADD_NEKTAR_TEST(PPF_R15000_ModifiedArnoldi_Shift)
    ADD_NEKTAR_TEST(ChanFlow_m8_Flowrate)
    ADD_NEKTAR_TEST(ChanFlow_3DH1D_FlowrateExplicit_MVM)
    ADD_NEKTAR_TEST(PhaseSampling)

    IF (NEKTAR_USE_SCOTCH)
        # This tests exhibits sensitivity to the choice of DirectStaticCond
        # vs. DirectMultiLevelStaticCond.
        ADD_NEKTAR_TEST(Pyr_channel_SVV)
    ENDIF()

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
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_FFT_ConOBC)
        ADD_NEKTAR_TEST(ChanFlow_3DH2D_FFT)
        ADD_NEKTAR_TEST(CylFlow_MovBody)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_MVM_FFTW_Consistency)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_6modes_FFTW_MixedDeal)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_P8_16modes_Mapping-explicit)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_P8_16modes_Mapping-implicit)
        ADD_NEKTAR_TEST(KovaFlow_3DH2D_xy_FFT)
        ADD_NEKTAR_TEST(KovaFlow_3DH2D_xz_FFT)
        ADD_NEKTAR_TEST(KovaFlow_m8_short_HOBC_3D1H)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_adaptive_16modes_FFTW_Mapping LENGTHY)
        ADD_NEKTAR_TEST(Couette_3DH2D_FFT)
    ENDIF (NEKTAR_USE_FFTW)

    IF (NEKTAR_USE_MPI)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Parallel_mode1)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Parallel_mode2)
        ADD_NEKTAR_TEST(ChanFlow_m3_par)
        ADD_NEKTAR_TEST(ChanFlow_m8_BodyForce_par LENGTHY)

        IF (NOT WIN32)
            ADD_NEKTAR_TEST(KovaFlow_expFromFile_par)
            ADD_NEKTAR_TEST(KovaFlow_varP_per_par)
            ADD_NEKTAR_TEST(Hex_channel_varP_par)
        ENDIF ()

        ADD_NEKTAR_TEST(Hex_channel_m3_par LENGTHY)
        ADD_NEKTAR_TEST(Pyr_channel_m6_par LENGTHY)
        ADD_NEKTAR_TEST(Hex_channel_m3_srhs_par LENGTHY)

        IF (NOT WIN32)
            ADD_NEKTAR_TEST(Tet_channel_m4_per_xxt_ml_par)
        ENDIF ()

        ADD_NEKTAR_TEST(Tet_channel_m8_par)
        ADD_NEKTAR_TEST(Tet_channel_m8_iter_ml_par LENGTHY)
        ADD_NEKTAR_TEST(bfs_tg_par LENGTHY)

        ADD_NEKTAR_TEST(ChanFlow_m8_Flowrate_par)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Flowrate_MVM)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Flowrate_MVM_hybrid)
        ADD_NEKTAR_TEST(Hex_channel_m3_Flowrate)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_FlowrateExplicit_MVM_par)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_FlowrateExplicit_MVM_par_hybrid)

        IF (NEKTAR_USE_PETSC)
            ADD_NEKTAR_TEST(Tet_channel_m8_petsc_sc_par LENGTHY)
        ENDIF (NEKTAR_USE_PETSC)
        IF (NEKTAR_USE_HDF5)
            ADD_NEKTAR_TEST(TriQuadChannelHdf5)
            ADD_NEKTAR_TEST(Tet_channel_m8_hdf)
            ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_MVM_hdf5)
            ADD_NEKTAR_TEST(CylFlow2D_FieldConvertFilter_Hdf5)
            IF (NEKTAR_USE_SCOTCH)
                ADD_NEKTAR_TEST(CubeAllElements_ChanFlow)
            ENDIF()
        ENDIF (NEKTAR_USE_HDF5)
    ENDIF (NEKTAR_USE_MPI)
ENDIF( NEKTAR_SOLVER_INCNAVIERSTOKES )
