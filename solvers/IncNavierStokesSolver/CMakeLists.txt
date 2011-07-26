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

IF (NEKTAR_USE_ARPACK)
    SET(IncNavierStokesSolverSource ${IncNavierStokesSolverSource}
        ../Auxiliary/DriverArnoldi.cpp)
ENDIF (NEKTAR_USE_ARPACK)

ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver 
			${IncNavierStokesSolverSource})
                        
