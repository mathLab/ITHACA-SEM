SET(EnsiteToXmlSources EnsiteToXml.cpp)
ADD_UTILITIES_EXECUTABLE(EnsiteToXml util-extra ${EnsiteToXmlSources})

IF (NEKTAR_USE_VTK)
    SET(VtkToFldSources VtkToFld.cpp)
    ADD_UTILITIES_EXECUTABLE(VtkToFld util-extra ${VtkToFldSources})
ENDIF (NEKTAR_USE_VTK)

