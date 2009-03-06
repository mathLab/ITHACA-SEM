SET(FHNSolverSource 
       ./FHNSolver.cpp 
       ./FHN.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(FHNSolver 
			${FHNSolverSource})
