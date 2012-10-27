SET(IncNavierStokesSolverSource 
       ./EquationSystems/CFLtester.cpp
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./AdvectionTerms/AdjointAdvection.cpp
       ./AdvectionTerms/AdvectionTerm.cpp
       ./AdvectionTerms/LinearisedAdvection.cpp
       ./AdvectionTerms/NavierStokesAdvection.cpp
       ./AdvectionTerms/SkewSymmetricAdvection.cpp
       ./AdvectionTerms/NoAdvection.cpp
       ./IncNavierStokesSolver.cpp 
       )

ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver solvers 
			${IncNavierStokesSolverSource})

ADD_TEST(NAME ChanFlowM3
         COMMAND Tester ${CMAKE_CURRENT_SOURCE_DIR}/Tests/ChanFlowM3.tst)