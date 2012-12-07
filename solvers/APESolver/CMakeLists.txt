SET(NEKTAR_SOLVER_APE ON CACHE BOOL
    "Build the APE solver.")

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
