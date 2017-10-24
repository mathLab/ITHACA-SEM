CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_CARDIAC_EP
    "Build the Cardiac electrophysiology solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_CARDIAC_EP OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_CARDIAC_EP )
    ADD_SOLVER_EXECUTABLE(CardiacEPSolver
        SOURCES CardiacEPSolver.cpp
        LIBRARY_SOURCES
        EquationSystems/Monodomain.cpp
        EquationSystems/Bidomain.cpp
        EquationSystems/BidomainRoth.cpp
        CellModels/CellModel.cpp
        CellModels/FitzhughNagumo.cpp
        CellModels/AlievPanfilov.cpp
        CellModels/CourtemancheRamirezNattel98.cpp
        CellModels/FentonKarma.cpp
        CellModels/LuoRudy91.cpp
        CellModels/Fox02.cpp
        CellModels/Winslow99.cpp
        CellModels/TenTusscher06.cpp
        Filters/FilterCheckpointCellModel.cpp
        Filters/FilterElectrogram.cpp
        Filters/FilterBenchmark.cpp
        Filters/FilterCellHistoryPoints.cpp
        Stimuli/Stimulus.cpp
        Stimuli/StimulusCircle.cpp
        Stimuli/StimulusRect.cpp
        Stimuli/StimulusPoint.cpp
        Stimuli/Protocol.cpp
        Stimuli/ProtocolSingle.cpp
        Stimuli/ProtocolS1.cpp
        Stimuli/ProtocolS1S2.cpp
        )

    ADD_SUBDIRECTORY(Utilities)
ENDIF( NEKTAR_SOLVER_CARDIAC_EP )
