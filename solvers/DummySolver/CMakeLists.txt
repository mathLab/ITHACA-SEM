CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_DUMMY
    "Build the Dummy solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_DUMMY )

    SET(DummySolverSource
       ./DummySolver.cpp
       ./EquationSystems/Dummy.cpp)

    ADD_SOLVER_EXECUTABLE(DummySolver SOURCES ${DummySolverSource})

#     ADD_NEKTAR_TEST(Dummy_3DCubeFile)
    ADD_NEKTAR_TEST(Dummy_3DCubeFileSingle)

    IF( NEKTAR_USE_CWIPI AND NEKTAR_USE_PETSC)
        ADD_NEKTAR_TEST(Dummy_3DCubeCwipi)
    ENDIF( NEKTAR_USE_CWIPI AND NEKTAR_USE_PETSC)

ENDIF( NEKTAR_SOLVER_DUMMY )
