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

    # This test requires MUMPS which is only built with FORTRAN support.
    IF( NEKTAR_USE_CWIPI AND NEKTAR_USE_PETSC AND CMAKE_Fortran_COMPILER)
        ADD_NEKTAR_TEST(Dummy_3DCubeCwipi)
    ENDIF()

ENDIF( NEKTAR_SOLVER_DUMMY )
