SET(FitzHughNagumoSolverSource 
       ./FitzHughNagumoSolver.cpp 
       ./FitzHughNagumo.cpp)

ADD_SOLVER_EXECUTABLE(FitzHughNagumoSolver solvers-extra 
			${FitzHughNagumoSolverSource})
