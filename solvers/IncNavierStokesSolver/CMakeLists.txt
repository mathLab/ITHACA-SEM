SET(IncNavierStokesSolverSource 
       ../Auxiliary/Driver.cpp
       ../Auxiliary/DriverStandard.cpp
       ../Auxiliary/DriverModifiedArnoldi.cpp
       ../Auxiliary/EquationSystem.cpp
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./EquationSystems/AdvectionTerm.cpp
       ./EquationSystems/NavierStokesAdvection.cpp
       ./EquationSystems/LinearisedAdvection.cpp
       ./EquationSystems/AdjointAdvection.cpp
       ./IncNavierStokesSolver.cpp 
       )

ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver 
			${IncNavierStokesSolverSource})
                        
