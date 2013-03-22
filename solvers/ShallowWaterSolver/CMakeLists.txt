CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_SHALLOW_WATER
    "Build the Shallow Water solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_SHALLOW_WATER )
    SET(ShallowWaterSolverSource 
       ./ShallowWaterSolver.cpp
       ./EquationSystems/ShallowWaterSystem.cpp
       ./EquationSystems/LinearSWE.cpp
       ./EquationSystems/NonlinearSWE.cpp)

    ADD_SOLVER_EXECUTABLE(ShallowWaterSolver solvers 
			${ShallowWaterSolverSource})

    ADD_NEKTAR_TEST(StandingWave_CG_P4)
    ADD_NEKTAR_TEST(StandingWave_CG_P8)
    ADD_NEKTAR_TEST(StandingWave_DG_P4)
    ADD_NEKTAR_TEST(StandingWave_DG_P8)
ENDIF( NEKTAR_SOLVER_SHALLOW_WATER )
