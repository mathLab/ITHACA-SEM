CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_FITZHUGH_NAGUMO
    "Build the Fitzhugh-Nagumo solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_FITZHUGH_NAGUMO )
    SET(FitzHughNagumoSolverSource 
       ./FitzHughNagumoSolver.cpp 
       ./FitzHughNagumo.cpp)

    ADD_SOLVER_EXECUTABLE(FitzHughNagumoSolver solvers-extra 
			${FitzHughNagumoSolverSource})
ENDIF( NEKTAR_SOLVER_FITZHUGH_NAGUMO )
