SET(IncNavierStokesSolverSource 
       ./IncNavierStokesSolver.cpp 
       ./IncNavierStokes.cpp
       ./VelocityCorrectionScheme.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver 
			${IncNavierStokesSolverSource})
