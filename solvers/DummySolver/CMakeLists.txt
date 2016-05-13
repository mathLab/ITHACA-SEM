CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_Dummy
    "Build the Dummy solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_Dummy AND NEKTAR_USE_CWIPI )

    SET(DummySolverSource
       ./DummySolver.cpp
       ./EquationSystems/Dummy.cpp)

    ADD_SOLVER_EXECUTABLE(DummySolver solvers ${DummySolverSource})

ENDIF( NEKTAR_SOLVER_Dummy AND NEKTAR_USE_CWIPI)
