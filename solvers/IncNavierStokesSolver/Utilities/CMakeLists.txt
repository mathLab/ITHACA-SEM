SET(CFLStepSources  CFLStep.cpp 
       ../EquationSystems/CoupledLinearNS.cpp
       ../EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ../EquationSystems/IncNavierStokes.cpp
       ../EquationSystems/VelocityCorrectionScheme.cpp
       ../AdvectionTerms/AdjointAdvection.cpp
       ../AdvectionTerms/AdvectionTerm.cpp
       ../AdvectionTerms/LinearisedAdvection.cpp
       ../AdvectionTerms/NavierStokesAdvection.cpp
       ../AdvectionTerms/SkewSymmetricAdvection.cpp
       ../AdvectionTerms/NoAdvection.cpp
       )

ADD_SOLVER_EXECUTABLE(CFLStep solvers ${CFLStepSources})

SET(AliasingSources  Aliasing.cpp 
       ../EquationSystems/CoupledLinearNS.cpp
       ../EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ../EquationSystems/IncNavierStokes.cpp
       ../EquationSystems/VelocityCorrectionScheme.cpp
       ../AdvectionTerms/AdjointAdvection.cpp
       ../AdvectionTerms/AdvectionTerm.cpp
       ../AdvectionTerms/LinearisedAdvection.cpp
       ../AdvectionTerms/NavierStokesAdvection.cpp
       ../AdvectionTerms/SkewSymmetricAdvection.cpp
       ../AdvectionTerms/NoAdvection.cpp
       )


ADD_SOLVER_EXECUTABLE(Aliasing solvers ${AliasingSources})
