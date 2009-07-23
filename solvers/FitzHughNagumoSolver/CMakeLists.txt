SET(FitzHughNagumoSolverSource 
       ./FitzHughNagumoSolver.cpp 
       ./FitzHughNagumo.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(FitzHughNagumoSolver 
			${FitzHughNagumoSolverSource})
