SET(VortexWaveInteractionSolverSource
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
        ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ../IncNavierStokesSolver/EquationSystems/CFLtester.cpp
       ../IncNavierStokesSolver/EquationSystems/CoupledLinearNS.cpp
       ../IncNavierStokesSolver/EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ../IncNavierStokesSolver/EquationSystems/IncNavierStokes.cpp
       ../IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp
       ../IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.cpp
       ../IncNavierStokesSolver/AdvectionTerms/NavierStokesAdvection.cpp
       ../IncNavierStokesSolver/AdvectionTerms/LinearisedAdvection.cpp
       ../IncNavierStokesSolver/AdvectionTerms/AdjointAdvection.cpp
       ./VortexWaveInteraction.cpp 
       ./VortexWaveInteractionSolver.cpp 
       )

SET(CalcL2ToLinfPressureSource
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ../IncNavierStokesSolver/EquationSystems/CFLtester.cpp
       ../IncNavierStokesSolver/EquationSystems/CoupledLinearNS.cpp
       ../IncNavierStokesSolver/EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ../IncNavierStokesSolver/EquationSystems/IncNavierStokes.cpp
       ../IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp
       ../IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.cpp
       ../IncNavierStokesSolver/AdvectionTerms/NavierStokesAdvection.cpp
       ../IncNavierStokesSolver/AdvectionTerms/LinearisedAdvection.cpp
       ../IncNavierStokesSolver/AdvectionTerms/AdjointAdvection.cpp
       ./CalcL2ToLinfPressure.cpp
       ./VortexWaveInteraction.h
       ./VortexWaveInteraction.cpp 
       )

SET(CalcVWIplusPressureSource
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ../IncNavierStokesSolver/EquationSystems/CFLtester.cpp
       ../IncNavierStokesSolver/EquationSystems/CoupledLinearNS.cpp
       ../IncNavierStokesSolver/EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ../IncNavierStokesSolver/EquationSystems/IncNavierStokes.cpp
       ../IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp
       ../IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.cpp
       ../IncNavierStokesSolver/AdvectionTerms/NavierStokesAdvection.cpp
       ../IncNavierStokesSolver/AdvectionTerms/LinearisedAdvection.cpp
       ../IncNavierStokesSolver/AdvectionTerms/AdjointAdvection.cpp
       ./CalcVWIplusPressure.cpp
       ./VortexWaveInteraction.cpp 
       )


SUBDIRS(utilities)


#IF (NEKTAR_USE_ARPACK)
#    SET(VortexWaveInteractionSolverSource ${VortexWaveInteractionSolverSource}
#        ../Auxiliary/DriverArpack.cpp)
#ENDIF (NEKTAR_USE_ARPACK)

ADD_SOLVER_EXECUTABLE(VortexWaveInteractionSolver solvers 
		     ${VortexWaveInteractionSolverSource})

ADD_SOLVER_EXECUTABLE(CalcL2ToLinfPressure solvers 
		     ${CalcL2ToLinfPressureSource})
                        
ADD_SOLVER_EXECUTABLE(CalcVWIplusPressure solvers 
		     ${CalcVWIplusPressureSource})
                        
