SET(LinkLibraries MultiRegions)

SET(PP_SOURCES ./Prepacing.cpp
        ../../CellModels/CellModel.cpp
        ../../CellModels/CourtemancheRamirezNattel98.cpp
        ../../CellModels/FentonKarma.cpp
        ../../CellModels/TenTusscher06.cpp
        ../../Stimuli/Stimulus.cpp
        ../../Stimuli/StimulusPoint.cpp
        ../../Stimuli/Protocol.cpp
        ../../Stimuli/ProtocolS1.cpp
        ../../Stimuli/ProtocolS1S2.cpp)
ADD_SOLVER_EXECUTABLE(PrePacing solvers-extra ${PP_SOURCES})

TARGET_LINK_LIBRARIES(PrePacing ${LinkLibraries})

ADD_NEKTAR_TEST(Courtemanche)
ADD_NEKTAR_TEST(CourtemancheAF)
ADD_NEKTAR_TEST(FentonKarma)
