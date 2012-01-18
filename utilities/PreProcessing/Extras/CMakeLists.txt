SET(EnsiteToXmlSources EnsiteToXml.cpp)
ADD_UTILITIES_EXECUTABLE(EnsiteToXml util-extra ${EnsiteToXmlSources})

IF (NEKTAR_USE_VTK)
    SET(VtkToFldSources VtkToFld.cpp)
    SET(VtkStripsToPolysSources VtkStripsToPolys.cpp)
    ADD_UTILITIES_EXECUTABLE(VtkToFld util-extra ${VtkToFldSources})
    ADD_UTILITIES_EXECUTABLE(VtkStripsToPolys util-extra ${VtkStripsToPolysSources})
ENDIF (NEKTAR_USE_VTK)

