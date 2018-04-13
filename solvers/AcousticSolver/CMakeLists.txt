CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_ACOUSTIC
    "Build the APE solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_ACOUSTIC OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_ACOUSTIC )
    SET(AcousticSolverSource 
       ./AcousticSolver.cpp
       ./EquationSystems/APE.cpp
       ./EquationSystems/LEE.cpp
       ./RiemannSolvers/AcousticSolver.cpp
       ./RiemannSolvers/LEESolver.cpp
       ./RiemannSolvers/APEUpwindSolver.cpp
       ./RiemannSolvers/APELaxFriedrichsSolver.cpp
       ./RiemannSolvers/LEEUpwindSolver.cpp
       ./RiemannSolvers/LEELaxFriedrichsSolver.cpp
)

    ADD_SOLVER_EXECUTABLE(AcousticSolver SOURCES ${AcousticSolverSource})

    ADD_NEKTAR_TEST(APE_1DPulseSource_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_1DPulseWall_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DChannel_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DPulseAdv_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DPulseInterp_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DPulseWall_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_2DVariableC_WeakDG_MODIFIED)
    ADD_NEKTAR_TEST(APE_3DPulse_WeakDG_MODIFIED LENGTHY)
    ADD_NEKTAR_TEST(APE_3DPulseWall_WeakDG_MODIFIED LENGTHY)

ENDIF( NEKTAR_SOLVER_ACOUSTIC )
