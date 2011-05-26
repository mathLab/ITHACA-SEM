SET(IncNavierStokesSolverSource 
       ../ADRSolver/EquationSystem.cpp
       ../ADRSolver/SessionReader.cpp
       ../Auxiliary/ADRBase.cpp
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./EquationSystems/AdvectionTerm.cpp
       ./EquationSystems/NavierStokesAdvection.cpp
       ./EquationSystems/LinearisedAdvection.cpp
       ./IncNavierStokesSolver.cpp 
       )

ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver 
			${IncNavierStokesSolverSource})
                        
