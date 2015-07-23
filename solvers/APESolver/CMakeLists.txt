CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_APE
    "Build the APE solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_APE )
    SET(APESolverSource 
       ./APESolver.cpp
       ./EquationSystems/APE.cpp
       ./RiemannSolvers/APESolver.cpp
       ./RiemannSolvers/UpwindSolver.cpp
       ./RiemannSolvers/LaxFriedrichsSolver.cpp)

    ADD_SOLVER_EXECUTABLE(APESolver solvers 
			${APESolverSource})

# TODO: We have a bug somewhere in the 1D boundary conditions. Therefore 1D
# problems are currently disabled. This should get fixed in the future. (See APE.cpp:64)
#    ADD_NEKTAR_TEST(APE_1DPulseSource_WeakDG_MODIFIED)
#    ADD_NEKTAR_TEST(APE_1DPulseWall_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DChannel_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DPulseAdv_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DPulseInterp_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DPulseWall_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DVariableC_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST_LENGTHY(APE_3DPulse_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST_LENGTHY(APE_3DPulseWall_WeakDG_MODIFIED)

ENDIF( NEKTAR_SOLVER_APE )
