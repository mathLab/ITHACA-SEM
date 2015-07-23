SET(CFLStepSources  CFLStep.cpp 
       ../EquationSystems/CoupledLinearNS.cpp
       ../EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ../EquationSystems/IncNavierStokes.cpp
       ../EquationSystems/VelocityCorrectionScheme.cpp
       ../EquationSystems/Extrapolate.cpp
       ../EquationSystems/StandardExtrapolate.cpp
       ../EquationSystems/SubSteppingExtrapolate.cpp
       ../AdvectionTerms/AdjointAdvection.cpp
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
       ../EquationSystems/Extrapolate.cpp
       ../EquationSystems/StandardExtrapolate.cpp
       ../EquationSystems/SubSteppingExtrapolate.cpp
       ../AdvectionTerms/AdjointAdvection.cpp
       ../AdvectionTerms/LinearisedAdvection.cpp
       ../AdvectionTerms/NavierStokesAdvection.cpp
       ../AdvectionTerms/SkewSymmetricAdvection.cpp
       ../AdvectionTerms/NoAdvection.cpp
       )


ADD_SOLVER_EXECUTABLE(Fld2DTo2D5            solvers Fld2DTo2D5.cpp)
ADD_SOLVER_EXECUTABLE(FldAddFalknerSkanBL   solvers FldAddFalknerSkanBL.cpp)
ADD_SOLVER_EXECUTABLE(AddModeTo2DFld        solvers AddModeTo2DFld.cpp)
ADD_SOLVER_EXECUTABLE(ExtractMeanModeFromHomo1DFld  
                                    solvers ExtractMeanModeFromHomo1DFld.cpp)
ADD_SOLVER_EXECUTABLE(Aliasing              solvers ${AliasingSources})
