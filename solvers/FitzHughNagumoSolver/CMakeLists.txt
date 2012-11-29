SET(NEKTAR_SOLVER_FITZHUGH_NAGUMO ON CACHE BOOL
    "Build the Fitzhugh-Nagumo solver.")

IF( NEKTAR_SOLVER_FITZHUGH_NAGUMO )
    SET(FitzHughNagumoSolverSource 
       ./FitzHughNagumoSolver.cpp 
       ./FitzHughNagumo.cpp)

    ADD_SOLVER_EXECUTABLE(FitzHughNagumoSolver solvers-extra 
			${FitzHughNagumoSolverSource})
ENDIF( NEKTAR_SOLVER_FITZHUGH_NAGUMO )
