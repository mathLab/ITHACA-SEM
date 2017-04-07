CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_VORTEXWAVE
    "Build the Vortex-Wave interaction solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_VORTEXWAVE )
    SET(VortexWaveInteractionSolverSource
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ./VortexWaveInteraction.cpp 
       ./VortexWaveInteractionSolver.cpp 
       $<TARGET_OBJECTS:IncNavierStokesSolver-obj>
       )

    SET(CalcL2ToLinfPressureSource
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ./CalcL2ToLinfPressure.cpp
       ./VortexWaveInteraction.h
       ./VortexWaveInteraction.cpp 
       $<TARGET_OBJECTS:IncNavierStokesSolver-obj>
       )

    SET(CalcVWIplusPressureSource
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ./CalcVWIplusPressure.cpp
       ./VortexWaveInteraction.cpp 
       $<TARGET_OBJECTS:IncNavierStokesSolver-obj>
       )

    SUBDIRS(utilities)

    ADD_SOLVER_EXECUTABLE(VortexWaveInteractionSolver SOURCES
	${VortexWaveInteractionSolverSource})
    ADD_SOLVER_UTILITY(CalcL2ToLinfPressure VortexWaveInteractionSolver SOURCES
	${CalcL2ToLinfPressureSource})
    ADD_SOLVER_UTILITY(CalcVWIplusPressure VortexWaveInteractionSolver SOURCES
	${CalcVWIplusPressureSource})

    FINALISE_CPACK_COMPONENT(vortexwaveinteraction-solver
        SUMMARY "Nektar++ vortex-wave interaction solver"
        DESCRIPTION "Nektar++ vortex-wave interaction solver")
ENDIF( NEKTAR_SOLVER_VORTEXWAVE )
