SET(IncNavierStokesSolverSource 
       ../Auxiliary/Driver.cpp
       ../Auxiliary/DriverStandard.cpp
       ../Auxiliary/DriverModifiedArnoldi.cpp
       ../Auxiliary/EquationSystem.cpp
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./AdvectionTerms/AdvectionTerm.cpp
       ./AdvectionTerms/NavierStokesAdvection.cpp
       ./AdvectionTerms/LinearisedAdvection.cpp
       ./AdvectionTerms/AdjointAdvection.cpp
       ./IncNavierStokesSolver.cpp 
       )

IF (NEKTAR_USE_ARPACK)
    SET(IncNavierStokesSolverSource ${IncNavierStokesSolverSource}
        ../Auxiliary/DriverArpack.cpp)
ENDIF (NEKTAR_USE_ARPACK)

ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver 
			${IncNavierStokesSolverSource})
                        
