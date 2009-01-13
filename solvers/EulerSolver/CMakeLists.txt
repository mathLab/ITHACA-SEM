SET(EulerSolverSource 
       ./EulerSolver.cpp 
       ./EulerEquations.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(EulerSolver 
			${EulerSolverSource})
