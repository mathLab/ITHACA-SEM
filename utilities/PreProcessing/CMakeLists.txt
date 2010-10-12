

#SET(MeshConvertSources MeshConvert.cpp  Gmsh.cpp TetGen.cpp)
#ADD_UTILITIES_EXECUTABLE(MeshConvert ${MeshConvertSources})

SET(EnsiteToXmlSources EnsiteToXml.cpp)
ADD_UTILITIES_EXECUTABLE(EnsiteToXml ${EnsiteToXmlSources})

SUBDIRS( MeshConvert )
