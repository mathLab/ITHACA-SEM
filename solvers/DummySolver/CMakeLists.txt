CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_DUMMY
    "Build the Dummy solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF
    "NEKTAR_USE_CWIPI" OFF)

IF( NEKTAR_SOLVER_DUMMY )

    SET(DummySolverSource
       ./DummySolver.cpp
       ./EquationSystems/Dummy.cpp)

    ADD_SOLVER_EXECUTABLE(DummySolver SOURCES ${DummySolverSource})

    ADD_NEKTAR_TEST(Dummy_3DCubeCwipi)

ENDIF( NEKTAR_SOLVER_DUMMY )
