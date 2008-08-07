
SET(LinkLibraries 

    optimized tinyxml debug tinyxml
)

SET(GmshSources Gmsh.cpp)
ADD_UTILITIES_EXECUTABLE(Gmsh ${GmshSources})
TARGET_LINK_LIBRARIES(Gmsh ${LinkLibraries})

SET(TINYXML_SOURCES
${TINYXML_SRC_DIR}/tinyxml.cpp
${TINYXML_SRC_DIR}/tinyxmlparser.cpp 
#${TINYXML_SRC_DIR}/xmltest.cpp
${TINYXML_SRC_DIR}/tinyxmlerror.cpp 
${TINYXML_SRC_DIR}/tinystr.cpp
)

SET(TINYXML_HEADERS
${TINYXML_INCLUDE_DIR}/tinystr.h
${TINYXML_INCLUDE_DIR}/tinyxml.h
)

ADD_UTILITIES_LIBRARY(tinyxml STATIC ${TINYXML_SOURCES} ${TINYXML_HEADERS})
# INSTALL(TARGETS tinyxml RUNTIME DESTINATION ${UTILITIES_BIN_DIR} OPTIONAL
#                         ARCHIVE DESTINATION ${UTILITIES_LIB_DIR} OPTIONAL
#                         LIBRARY DESTINATION ${UTILITIES_LIB_DIR} OPTIONAL)
