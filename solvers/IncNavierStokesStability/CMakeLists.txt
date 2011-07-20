IF(NEKTAR_USE_ARPACK)
SET(IncNavierStokesStabilitySource 
       ../Auxiliary/EquationSystem.cpp
       ../Auxiliary/ADRBase.cpp
       ../IncNavierStokesSolver/EquationSystems/CoupledLinearNS.cpp
       ../IncNavierStokesSolver/EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ../IncNavierStokesSolver/EquationSystems/IncNavierStokes.cpp
       ../IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp
       ../IncNavierStokesSolver//EquationSystems/AdvectionTerm.cpp
       ../IncNavierStokesSolver/EquationSystems/NavierStokesAdvection.cpp
       ../IncNavierStokesSolver/EquationSystems/LinearisedAdvection.cpp
       ../IncNavierStokesSolver/EquationSystems/AdjointAdvection.cpp
       ./IncNavierStokesStability.cpp 
       )

ADD_SOLVER_EXECUTABLE(IncNavierStokesStability
			${IncNavierStokesStabilitySource})

ENDIF(NEKTAR_USE_ARPACK)
