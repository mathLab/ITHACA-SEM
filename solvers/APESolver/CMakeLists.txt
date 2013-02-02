CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_APE
    "Build the APE solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_APE )
    SET(APESolverSource 
       ./APESolver.cpp
       ./EquationSystems/APESystem.cpp
       ./EquationSystems/APE.cpp)

    ADD_SOLVER_EXECUTABLE(APESolver solvers 
			${APESolverSource})

    ADD_NEKTAR_TEST(Channel)
    ADD_NEKTAR_TEST(Pulse)
ENDIF( NEKTAR_SOLVER_APE )
