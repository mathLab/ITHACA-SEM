SET(FitzHughNagumoSolverSource 
       ./FitzHughNagumoSolver.cpp 
       ./FitzHughNagumo.cpp
       ../Auxiliary/EquationSystem.cpp)

ADD_SOLVER_EXECUTABLE(FitzHughNagumoSolver 
			${FitzHughNagumoSolverSource})
