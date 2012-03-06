SET(VortexWaveInteractionSolverSource
#       ../Auxiliary/Driver.cpp
#       ../Auxiliary/DriverStandard.cpp
#       ../Auxiliary/DriverArnoldi.cpp
#       ../Auxiliary/DriverModifiedArnoldi.cpp
#       ../Auxiliary/EquationSystem.cpp
#       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
#       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
#       ../IncNavierStokesSolver/EquationSystems/CoupledLinearNS.cpp
#       ../IncNavierStokesSolver/EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
#       ../IncNavierStokesSolver/EquationSystems/IncNavierStokes.cpp
#       ../IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp
#       ../IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.cpp
#       ../IncNavierStokesSolver/AdvectionTerms/NavierStokesAdvection.cpp
#       ../IncNavierStokesSolver/AdvectionTerms/LinearisedAdvection.cpp
#       ../IncNavierStokesSolver/AdvectionTerms/AdjointAdvection.cpp
#       ./VortexWaveInteraction.cpp 
#       ./VortexWaveInteractionSolver.cpp 
       )

IF (NEKTAR_USE_ARPACK)
    SET(VortexWaveInteractionSolverSource ${VortexWaveInteractionSolverSource}
        ../Auxiliary/DriverArpack.cpp)
ENDIF (NEKTAR_USE_ARPACK)

#ADD_SOLVER_EXECUTABLE(VortexWaveInteractionSolver solvers 
#		     ${VortexWaveInteractionSolverSource})
                        
