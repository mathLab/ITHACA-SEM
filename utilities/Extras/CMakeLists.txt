ADD_UTILITIES_EXECUTABLE(CheckXmlFile COMPONENT utilities-extra SOURCES CheckXmlFile.cpp)
ADD_UTILITIES_EXECUTABLE(StreamFunction2D COMPONENT utilities-extra SOURCES StreamFunction2D.cpp)

IF (NEKTAR_USE_VTK)
    ADD_UTILITIES_EXECUTABLE(VtkToFld COMPONENT utilities-extra SOURCES VtkToFld.cpp)
    ADD_UTILITIES_EXECUTABLE(VtkStripsToPolys COMPONENT utilities-extra SOURCES VtkStripsToPolys.cpp)

    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(VtkToFld LINK_PUBLIC vtkCommon vtkIO vtkGraphics)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys LINK_PUBLIC vtkCommon vtkIO vtkGraphics)
    ELSEIF(VTK_MAJOR_VERSION GREATER 6)
        TARGET_LINK_LIBRARIES(VtkToFld LINK_PUBLIC vtkFiltersCore vtkFiltersParallelDIY2 vtkCommonCore vtkIOLegacy
                                               vtkFiltersGeometry)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys LINK_PUBLIC vtkCommonCore vtkIOLegacy)
    ELSE ()
        TARGET_LINK_LIBRARIES(VtkToFld LINK_PUBLIC vtkCommonCore vtkIOLegacy
                                                   vtkFiltersGeometry)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys LINK_PUBLIC vtkCommonCore vtkIOLegacy)
    ENDIF ()

    ADD_UTILITIES_EXECUTABLE(VtkToPng COMPONENT utilities-extra SOURCES VtkToPng.cpp)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(VtkToPng LINK_PUBLIC vtkCommon vtksys vtkViews vtkWidgets
                                                   vtkRendering vtkIO)
    ELSE ()
        TARGET_LINK_LIBRARIES(VtkToPng LINK_PUBLIC ${VTK_LIBRARIES})
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)
