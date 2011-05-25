SET(IncNavierStokesSolverSource 
       ../ADRSolver/EquationSystem.cpp
       ../ADRSolver/SessionReader.cpp
       ../Auxiliary/ADRBase.cpp
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./IncNavierStokesSolver.cpp 
       )

ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver 
			${IncNavierStokesSolverSource})
                        
