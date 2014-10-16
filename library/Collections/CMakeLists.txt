SET(COLLECTIONS_SOURCES
  Collection.cpp
  Operator.cpp
  HexOperators.cpp
  IProduct.cpp
  PrismOperators.cpp
  QuadOperators.cpp
  SegOperators.cpp
  TetOperators.cpp
  TriOperators.cpp
  )

SET(SOLVER_UTILS_HEADERS
  Collection.h
  IProduct.h
  Operator.h
)

#ADD_DEFINITIONS(-DSOLVER_UTILS_EXPORTS)

ADD_NEKTAR_LIBRARY(Collections lib ${NEKTAR_LIBRARY_TYPE} ${COLLECTIONS_SOURCES} ${COLLECTIONS_HEADERS})
TARGET_LINK_LIBRARIES(Collections LINK_PUBLIC LocalRegions)

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/Collections COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
