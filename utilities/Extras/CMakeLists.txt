ADD_UTILITIES_EXECUTABLE(CheckXmlFile COMPONENT utilities-extra SOURCES CheckXmlFile.cpp)
ADD_UTILITIES_EXECUTABLE(StreamFunction2D COMPONENT utilities-extra SOURCES StreamFunction2D.cpp)

IF (NEKTAR_USE_VTK)
    ADD_UTILITIES_EXECUTABLE(VtkToFld COMPONENT utilities-extra SOURCES VtkToFld.cpp)
    ADD_UTILITIES_EXECUTABLE(VtkStripsToPolys COMPONENT utilities-extra SOURCES VtkStripsToPolys.cpp)

    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(VtkToFld vtkCommon vtkIO vtkGraphics)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys vtkCommon vtkIO vtkGraphics)
    ELSE ()
        TARGET_LINK_LIBRARIES(VtkToFld vtkCommonCore vtkIOLegacy
                                       vtkFiltersGeometry)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys vtkCommonCore vtkIOLegacy)
    ENDIF ()
    
    ADD_UTILITIES_EXECUTABLE(VtkToPng COMPONENT utilities-extra SOURCES VtkToPng.cpp)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(VtkToPng vtkCommon vtksys vtkViews vtkWidgets 
                                       vtkRendering vtkIO)
    ELSE ()
        TARGET_LINK_LIBRARIES(VtkToPng ${VTK_LIBRARIES})
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)

