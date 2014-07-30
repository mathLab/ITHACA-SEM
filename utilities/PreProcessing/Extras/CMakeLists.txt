SET(EnsiteToXmlSources EnsiteToXml.cpp)
ADD_UTILITIES_EXECUTABLE(EnsiteToXml util-extra ${EnsiteToXmlSources})

SET(CheckXmlFileSources CheckXmlFile.cpp)
ADD_UTILITIES_EXECUTABLE(CheckXmlFile util-extra ${CheckXmlFileSources})

IF (NEKTAR_USE_VTK)
    SET(VtkToFldSources VtkToFld.cpp)
    SET(VtkStripsToPolysSources VtkStripsToPolys.cpp)
    ADD_UTILITIES_EXECUTABLE(VtkToFld util-extra ${VtkToFldSources})
    ADD_UTILITIES_EXECUTABLE(VtkStripsToPolys util-extra ${VtkStripsToPolysSources})
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(VtkToFld vtkCommon vtkIO vtkGraphics)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys vtkCommon vtkIO vtkGraphics)
    ELSE ()
        TARGET_LINK_LIBRARIES(VtkToFld vtkCommonCore vtkIOLegacy
                                       vtkFiltersGeometry)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys vtkCommonCore vtkIOLegacy)
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)

