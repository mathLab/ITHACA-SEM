CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_CARDIAC_EP
    "Build the Cardiac electrophysiology solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_CARDIAC_EP )
    SET(CardiacEPSolverSource 
       ./CardiacEPSolver.cpp
       ./EquationSystems/Monodomain.cpp
       ./EquationSystems/Bidomain.cpp
       ./EquationSystems/BidomainRoth.cpp
       ./CellModels/CellModel.cpp
       ./CellModels/FitzhughNagumo.cpp
       ./CellModels/AlievPanfilov.cpp
       ./CellModels/CourtemancheRamirezNattel98.cpp
       ./CellModels/FentonKarma.cpp
       ./CellModels/LuoRudy91.cpp
       ./CellModels/Fox02.cpp
       ./CellModels/Winslow99.cpp
       ./CellModels/TenTusscher06.cpp
       ./Filters/FilterCheckpointCellModel.cpp
       ./Filters/FilterElectrogram.cpp
       ./Filters/FilterBenchmark.cpp
       ./Stimuli/Stimulus.cpp
       ./Stimuli/StimulusCircle.cpp
       ./Stimuli/StimulusRect.cpp
       ./Stimuli/StimulusPoint.cpp
       ./Stimuli/Protocol.cpp
       ./Stimuli/ProtocolSingle.cpp
       ./Stimuli/ProtocolS1.cpp
       ./Stimuli/ProtocolS1S2.cpp
    )

    ADD_SOLVER_EXECUTABLE(CardiacEPSolver solvers-extra 
			${CardiacEPSolverSource})

    ADD_SUBDIRECTORY(Utilities)
ENDIF( NEKTAR_SOLVER_CARDIAC_EP )
