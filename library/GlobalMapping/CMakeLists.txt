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

ADD_NEKTAR_LIBRARY(GlobalMapping
    SOURCES ${GLOBAL_MAPPING_SOURCES}
    HEADERS ${GLOBAL_MAPPING_HEADERS}
    DEPENDS MultiRegions
    SUMMARY "Nektar++ GlobalMapping library"
    DESCRIPTION "This library provides deformation routines for Nektar++ simulations.")
