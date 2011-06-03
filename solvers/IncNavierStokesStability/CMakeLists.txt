IF(NEKTAR_USE_ARPACK)
SET(IncNavierStokesStabilitySource 
       ../ADRSolver/EquationSystem.cpp
       ../Auxiliary/ADRBase.cpp
       ../IncNavierStokesSolver/EquationSystems/CoupledLinearNS.cpp
       ../IncNavierStokesSolver/EquationSystems/IncNavierStokes.cpp
       ../IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp
       ../IncNavierStokesSolver//EquationSystems/AdvectionTerm.cpp
       ../IncNavierStokesSolver/EquationSystems/NavierStokesAdvection.cpp
       ../IncNavierStokesSolver/EquationSystems/LinearisedAdvection.cpp
       ./IncNavierStokesStability.cpp 
       )

ADD_SOLVER_EXECUTABLE(IncNavierStokesStability
			${IncNavierStokesStabilitySource})

ENDIF(NEKTAR_USE_ARPACK)
