SET(GLOBAL_MAPPING_SOURCES
  Mapping.cpp
  MappingXofXZ.cpp
  MappingXofZ.cpp
  MappingXYofZ.cpp
  MappingXYofT.cpp
  MappingXYofXY.cpp
  MappingGeneralFixed.cpp
  MappingIdentity.cpp
)

SET(GLOBAL_MAPPING_HEADERS
  Mapping.h
  MappingXofXZ.h
  MappingXofZ.h
  MappingXYofZ.h
  MappingXYofT.h
  MappingXYofXY.h
  MappingGeneralFixed.h
  MappingIdentity.h
)

ADD_DEFINITIONS(-DGLOBAL_MAPPING_EXPORTS)

ADD_NEKTAR_LIBRARY(GlobalMapping lib ${NEKTAR_LIBRARY_TYPE} ${GLOBAL_MAPPING_SOURCES} ${GLOBAL_MAPPING_HEADERS})
TARGET_LINK_LIBRARIES(GlobalMapping MultiRegions)

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/GlobalMapping COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
