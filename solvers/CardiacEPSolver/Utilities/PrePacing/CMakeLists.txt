SET(LinkLibraries
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g 
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug                             ${Boost_THREAD_LIBRARY_DEBUG}
    optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug                          ${Boost_IOSTREAMS_LIBRARY_DEBUG}
    optimized ${ZLIB_LIBRARY_RELEASE} debug ${ZLIB_LIBRARY_DEBUG}
    optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
    optimized sbtk    debug sbtk
)

SET(PP_SOURCES ./Prepacing.cpp
        ../../CellModels/CellModel.cpp
        ../../CellModels/CourtemancheRamirezNattel98.cpp
        ../../Stimuli/Stimulus.cpp
        ../../Stimuli/StimulusPoint.cpp
        ../../Stimuli/Protocol.cpp
        ../../Stimuli/ProtocolS1S2.cpp)
ADD_SOLVER_EXECUTABLE(PrePacing solvers-extra ${PP_SOURCES})
TARGET_LINK_LIBRARIES(PrePacing ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(PrePacing)
