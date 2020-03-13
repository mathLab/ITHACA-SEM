CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_ELASTICITY
  "Build the linear elasticity solver." ON
  "NEKTAR_BUILD_SOLVERS" OFF)

IF (NEKTAR_SOLVER_ELASTICITY)
    SET(LinearElasticSolverSource
      ./LinearElasticSolver.cpp 
      ./EquationSystems/LinearElasticSystem.cpp
      ./EquationSystems/IterativeElasticSystem.cpp
      ./EquationSystems/CoupledAssemblyMap.cpp)

    ADD_SOLVER_EXECUTABLE(LinearElasticSolver solvers 
	${LinearElasticSolverSource})

    ADD_NEKTAR_TEST(L-domain)
    IF (NEKTAR_USE_MPI)
        ADD_NEKTAR_TEST(L-domain-par)
    ENDIF (NEKTAR_USE_MPI)
ENDIF (NEKTAR_SOLVER_ELASTICITY)
