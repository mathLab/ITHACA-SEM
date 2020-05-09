ADD_NEKTAR_EXECUTABLE(Helmholtz1D
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz2D
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3D
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz3D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3DHomo1D
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz3DHomo1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3DHomo2D
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz3DHomo2D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz1D
    COMPONENT demos DEPENDS MultiRegions SOURCES HDGHelmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz2D
    COMPONENT demos DEPENDS MultiRegions SOURCES HDGHelmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz3D
    COMPONENT demos DEPENDS MultiRegions SOURCES HDGHelmholtz3D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz3DHomo1D
    COMPONENT demos DEPENDS MultiRegions SOURCES HDGHelmholtz3DHomo1D.cpp)
ADD_NEKTAR_EXECUTABLE(PostProcHDG2D
    COMPONENT demos DEPENDS MultiRegions SOURCES PostProcHDG2D.cpp)
ADD_NEKTAR_EXECUTABLE(PostProcHDG3D
    COMPONENT demos DEPENDS MultiRegions SOURCES PostProcHDG3D.cpp)
ADD_NEKTAR_EXECUTABLE(Deriv3DHomo1D_SingleMode
    COMPONENT demos DEPENDS MultiRegions SOURCES Deriv3DHomo1D_SingleMode.cpp)
ADD_NEKTAR_EXECUTABLE(Deriv3DHomo1D
    COMPONENT demos DEPENDS MultiRegions SOURCES Deriv3DHomo1D.cpp)
ADD_NEKTAR_EXECUTABLE(Deriv3DHomo2D
    COMPONENT demos DEPENDS MultiRegions SOURCES Deriv3DHomo2D.cpp)
ADD_NEKTAR_EXECUTABLE(SteadyAdvectionDiffusionReaction2D
    COMPONENT demos DEPENDS MultiRegions SOURCES SteadyAdvectionDiffusionReaction2D.cpp)

# Add ExtraDemos subdirectory.
IF (IS_DIRECTORY ExtraDemos)
    ADD_SUBDIRECTORY(ExtraDemos)
ENDIF()

ADD_NEKTAR_TEST(Helmholtz1D_CG_P8)
ADD_NEKTAR_TEST(Helmholtz1D_CG_P8_periodic)
ADD_NEKTAR_TEST(Helmholtz1D_CG_P8_RBC)
ADD_NEKTAR_TEST(Helmholtz1D_HDG_P8)
ADD_NEKTAR_TEST(Helmholtz1D_HDG_P8_RBC)

ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Nodes)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_sc)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_full)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_iter_sc)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P9_Modes_varcoeff)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_curved_quad)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_curved_tri)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_PreconBlock)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_PreconDiagonal)
ADD_NEKTAR_TEST(Helmholtz2D_HDG_P7_Modes)
ADD_NEKTAR_TEST(Helmholtz2D_HDG_P7_Modes_AllBCs)
ADD_NEKTAR_TEST(Helmholtz2D_CG_varP_Modes)

ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex LENGTHY)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_Heterogeneous)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_Heterogeneous_iterLE)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_iter_sc_cont)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Tet)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Tet_VarP)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Tet_VarP_iterLE)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism_VarP)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism_VarP_iterLE)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism_Deformed)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Pyr)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Pyr_VarP)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Pyr_VarP_iterLE)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Pyr_Deformed)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Homo1D)
ADD_NEKTAR_TEST(Helmholtz3D_HDG_Homo1D LENGTHY)
ADD_NEKTAR_TEST(Helmholtz3D_HDG_Prism)
ADD_NEKTAR_TEST(Helmholtz3D_HDG_Hex_AllBCs LENGTHY)
ADD_NEKTAR_TEST(Helmholtz3D_HDG_Tet)

ADD_NEKTAR_TEST(LinearAdvDiffReact2D_P7_Modes)
ADD_NEKTAR_TEST(Deriv3D_Homo1D)
ADD_NEKTAR_TEST(Deriv3D_Homo2D)

# Tests exclusively for multi-level static condensation (which will not run
# using Windows).
IF (NEKTAR_USE_SCOTCH)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_mlsc)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_iter_ml)
    ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism_iter_ml)
    ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_iter_ml)
ENDIF()

# PETSc solver tests
IF (NEKTAR_USE_PETSC)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_petsc_full)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_petsc_ml)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_petsc_sc)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_petsc_sc_shell)
ENDIF (NEKTAR_USE_PETSC)

IF (NEKTAR_USE_MPI)
    ADD_NEKTAR_TEST(CubeAllElements LENGTHY)
    ADD_NEKTAR_TEST(Helmholtz1D_CG_P8_periodic_par3)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_xxt_full)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_xxt_sc)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_iter_full_par2)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_iter_sc_par2)
    # These tests have been disabled on WIN32 because there's no support
    # for MLSC on Windows builds.
    IF (NOT WIN32)
        ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_iter_ml_par3)
        ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_xxt_ml)
        ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_iter_ml_par3)
    ENDIF ()
    ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism_iter_ml_par5)
    ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_xxt_sc_par3 LENGTHY)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P14_xxt_per)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_varP_Modes_par)

    # TODO: This test fails due to a bug with Dirichlet bnd conditions.
    #       To be resolved in a separate branch.
    # ADD_NEKTAR_TEST(Helmholtz3D_CG_Tet_iter_global_par3)

    ADD_NEKTAR_TEST(Helmholtz2D_HDG_P7_Modes_AllBCs_par2 LENGTHY)
    ADD_NEKTAR_TEST(Helmholtz3D_HDG_Prism_par2)
    ADD_NEKTAR_TEST(Helmholtz3D_HDG_Hex_AllBCs_par2)

    IF (NEKTAR_USE_METIS)
        IF (NOT WIN32)
            ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_iter_ml_metis_par3)
        ENDIF ()
        ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism_iter_ml_metis_par3)
        ADD_NEKTAR_TEST(Helmholtz3D_HDG_Hex_AllBCs_metis_par2)
    ENDIF (NEKTAR_USE_METIS)

    IF (NEKTAR_USE_PETSC)
        ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_petsc_sc_par3)
        ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_petsc_sc_par3)
        ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_petsc_sc_shell_par3)
    ENDIF (NEKTAR_USE_PETSC)
ENDIF (NEKTAR_USE_MPI)
