
SET(ExtractCriticalLayerSource
       ./ExtractCriticalLayer.cpp
       )


SET(MoveMeshToCriticalLayerSource
       ./MoveMeshToCriticalLayer.cpp
       )


SET(ExpandMeshByRotationSource
       ./ExpandMeshByRotation.cpp
       )

SET(MeshMoveSource
       ./MeshMove.cpp
       )
ADD_SOLVER_EXECUTABLE(ExtractCriticalLayer solvers 
		     ${ExtractCriticalLayerSource})

ADD_SOLVER_EXECUTABLE(MoveMeshToCriticalLayer solvers 
		     ${MoveMeshToCriticalLayerSource})

ADD_SOLVER_EXECUTABLE(ExpandMeshByRotation solvers 
		     ${ExpandMeshByRotationSource})

ADD_SOLVER_EXECUTABLE(MeshMove solvers 
		     ${MeshMoveSource})
                        
