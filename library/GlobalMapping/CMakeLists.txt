SET(GLOBAL_MAPPING_SOURCES
  Deform.cpp
  Mapping.cpp
  MappingXofXZ.cpp
  MappingXofZ.cpp
  MappingXYofZ.cpp
  MappingXYofXY.cpp
  MappingGeneral.cpp
  MappingTranslation.cpp
)

SET(GLOBAL_MAPPING_HEADERS
  Deform.h
  Mapping.h
  MappingXofXZ.h
  MappingXofZ.h
  MappingXYofZ.h
  MappingXYofXY.h
  MappingGeneral.h
  MappingTranslation.h
)

ADD_DEFINITIONS(-DGLOBAL_MAPPING_EXPORTS)

ADD_NEKTAR_LIBRARY(GlobalMapping lib ${NEKTAR_LIBRARY_TYPE} ${GLOBAL_MAPPING_SOURCES} ${GLOBAL_MAPPING_HEADERS})
TARGET_LINK_LIBRARIES(GlobalMapping MultiRegions)

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/GlobalMapping COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
