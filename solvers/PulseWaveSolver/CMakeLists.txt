CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_PULSEWAVE
    "Build the Pulse-wave solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_PULSEWAVE OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_PULSEWAVE )
    ADD_SOLVER_EXECUTABLE(PulseWaveSolver
        SOURCES PulseWaveSolver.cpp
        LIBRARY_SOURCES
        EquationSystems/PulseWaveSystem.h
        EquationSystems/PulseWaveSystem.cpp
        EquationSystems/PulseWavePropagation.h
        EquationSystems/PulseWavePropagation.cpp
        EquationSystems/PulseWaveBoundary.h
        EquationSystems/PulseWaveBoundary.cpp
        EquationSystems/PulseWavePressureArea.h
        EquationSystems/PulseWavePressureArea.cpp
        EquationSystems/TimeDependentInflow.h
        EquationSystems/TimeDependentInflow.cpp
        EquationSystems/QInflow.h
        EquationSystems/QInflow.cpp
        EquationSystems/AInflow.h
        EquationSystems/AInflow.cpp
        EquationSystems/UInflow.h
        EquationSystems/UInflow.cpp
        EquationSystems/TerminalOutflow.h
        EquationSystems/TerminalOutflow.cpp
        EquationSystems/ROutflow.h
        EquationSystems/ROutflow.cpp
        EquationSystems/RCROutflow.h
        EquationSystems/RCROutflow.cpp
        EquationSystems/UndefinedInOutflow.h
        EquationSystems/UndefinedInOutflow.cpp
        EquationSystems/LymphaticPressureArea.h
        EquationSystems/LymphaticPressureArea.cpp
        EquationSystems/ArterialPressureArea.h
        EquationSystems/ArterialPressureArea.cpp
        RiemannSolvers/UpwindPulseSolver.cpp
        )

    SUBDIRS(Utilities)

    ADD_NEKTAR_TEST(55_Artery_Network)
    ADD_NEKTAR_TEST(Bifurcation)
    ADD_NEKTAR_TEST(Junction)
    ADD_NEKTAR_TEST(TwoBifurcations)
    ADD_NEKTAR_TEST(Merging)
    ADD_NEKTAR_TEST(VariableAreaTest)
    ADD_NEKTAR_TEST(VariableMatPropTest)
    ADD_NEKTAR_TEST(Q_inflow)
    ADD_NEKTAR_TEST(A_inflow)
    ADD_NEKTAR_TEST(U_inflow)
    ADD_NEKTAR_TEST(RCR_boundary)
ENDIF( NEKTAR_SOLVER_PULSEWAVE )
