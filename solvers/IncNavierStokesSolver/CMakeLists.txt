SET(IncNavierStokesSolverSource 
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./AdvectionTerms/AdvectionTerm.cpp
       ./AdvectionTerms/NavierStokesAdvection.cpp
       ./AdvectionTerms/LinearisedAdvection.cpp
       ./AdvectionTerms/AdjointAdvection.cpp
       ./AdvectionTerms/SkewSymmetricAdvection.cpp
       ./IncNavierStokesSolver.cpp 
       )

ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver solvers 
			${IncNavierStokesSolverSource})
                        
