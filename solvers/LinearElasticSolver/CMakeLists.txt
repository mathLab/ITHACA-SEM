CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_ELASTICITY
  "Build the linear elasticity solver." ON
  "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_ELASTICITY OFF CACHE INTERNAL "")
ENDIF()

IF (NEKTAR_SOLVER_ELASTICITY)
    SET(LinearElasticSolverSource
      ./LinearElasticSolver.cpp 
      ./EquationSystems/LinearElasticSystem.cpp
      ./EquationSystems/IterativeElasticSystem.cpp
      ./EquationSystems/CoupledAssemblyMap.cpp)

    ADD_SOLVER_EXECUTABLE(LinearElasticSolver SOURCES ${LinearElasticSolverSource})

    ADD_NEKTAR_TEST(L-domain)
    IF (NEKTAR_USE_MPI)
        IF (WIN32)
            TARGET_LINK_LIBRARIES(LinearElasticSolver LINK_PRIVATE ${GSMPI_LIBRARY} ${XXT_LIBRARY})
        ENDIF (WIN32)
        ADD_NEKTAR_TEST(L-domain-par)
    ENDIF (NEKTAR_USE_MPI)
ENDIF (NEKTAR_SOLVER_ELASTICITY)
