ADD_SOLVER_UTILITY(ExtractCriticalLayer VortexWaveInteraction
    SOURCES ExtractCriticalLayerFunctions.cpp ExtractCriticalLayer.cpp)
ADD_SOLVER_UTILITY(MoveMeshToCriticalLayer VortexWaveInteraction
    SOURCES ExtractCriticalLayerFunctions.cpp MoveMeshToCriticalLayer.cpp)
ADD_SOLVER_UTILITY(ExpandMeshByRotation VortexWaveInteraction
    SOURCES ExpandMeshByRotation.cpp)
ADD_SOLVER_UTILITY(MeshMove VortexWaveInteraction
    SOURCES MeshMove.cpp)
