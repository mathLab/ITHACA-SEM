SET(COLLECTIONS_SOURCES
  Operator.cpp
  QuadOperators.cpp
  )

SET(SOLVER_UTILS_HEADERS
  Collection.h
  Operator.h
)

#ADD_DEFINITIONS(-DSOLVER_UTILS_EXPORTS)

ADD_NEKTAR_LIBRARY(Collections lib ${NEKTAR_LIBRARY_TYPE} ${COLLECTIONS_SOURCES} ${COLLECTIONS_HEADERS})
TARGET_LINK_LIBRARIES(Collections
  optimized LibUtilities   debug LibUtilities-g
  optimized SpatialDomains debug SpatialDomains-g
  optimized StdRegions     debug StdRegions-g
  optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
  ${Boost_THREAD_LIBRARY}
  ${Boost_IOSTREAMS_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${ZLIB_LIBRARY}
  )
SET_LAPACK_LINK_LIBRARIES(Collections)

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/Collections COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
