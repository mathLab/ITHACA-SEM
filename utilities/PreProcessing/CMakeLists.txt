
SET(LinkLibraries 

    optimized tinyxml debug tinyxml
)


SET(PreProcessingHeaders
     )

SET(PreProcessingSources
     Gmsh.cpp
     TetGen.cpp
     )
 ADD_UTILITIES_LIBRARY(preprocessing STATIC ${PreProcessingSources} ${PreProcessingHeaders})
 TARGET_LINK_LIBRARIES(preprocessing ${LinkLibraries})


# SET(GmshSources Gmsh.cpp)
# ADD_UTILITIES_EXECUTABLE(Gmsh ${GmshSources})
# TARGET_LINK_LIBRARIES(Gmsh ${LinkLibraries})
#
# SET(TetGenSources TetGen.cpp)
# ADD_UTILITIES_EXECUTABLE(TetGen ${TetGenSources})
# TARGET_LINK_LIBRARIES(TetGen ${LinkLibraries})


SET(MeshConvertSources MeshConvert.cpp)
ADD_UTILITIES_EXECUTABLE(MeshConvert ${MeshConvertSources})
TARGET_LINK_LIBRARIES(MeshConvert ${LinkLibraries})


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
