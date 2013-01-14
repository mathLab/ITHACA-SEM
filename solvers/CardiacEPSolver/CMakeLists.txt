CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_CARDIAC_EP
    "Build the Cardiac electrophysiology solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_CARDIAC_EP )
    SET(CardiacEPSolverSource 
       ./CardiacEPSolver.cpp
       ./EquationSystems/Monodomain.cpp
       ./EquationSystems/Bidomain.cpp
       ./CellModels/CellModel.cpp
       ./CellModels/FitzhughNagumo.cpp
       ./CellModels/AlievPanfilov.cpp
       ./CellModels/CourtemancheRamirezNattel98.cpp
       ./CellModels/LuoRudy91.cpp
       ./CellModels/Fox02.cpp
       ./CellModels/Winslow99.cpp
       ./CellModels/TenTusscher06Epi.cpp
       ./CellModels/TenTusscher06M.cpp
       ./CellModels/TenTusscher06Endo.cpp
       ./Filters/FilterCheckpointCellModel.cpp
       ./Stimuli/Stimulus.cpp
       ./Stimuli/StimulusCircle.cpp
       ./Stimuli/StimulusRect.cpp
       ./Stimuli/Protocol.cpp
       ./Stimuli/ProtocolSingle.cpp
       ./Stimuli/ProtocolS1.cpp
       ./Stimuli/ProtocolS1S2.cpp
    )

    ADD_SOLVER_EXECUTABLE(CardiacEPSolver solvers-extra 
			${CardiacEPSolverSource})
ENDIF( NEKTAR_SOLVER_CARDIAC_EP )
