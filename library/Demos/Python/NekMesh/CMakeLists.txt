ADD_NEKPY_EXECUTABLE(StructuredGrid.py StructuredGrid.py)

IF (NEKTAR_BUILD_TESTS)
    ADD_NEKPY_TEST(structured_grid)
ENDIF()

IF (NEKTAR_USE_MESHGEN)
    ADD_NEKPY_EXECUTABLE(LoadCAD.py LoadCAD.py)

    IF (NEKTAR_BUILD_TESTS)
        ADD_NEKPY_TEST(load_cad)
    ENDIF()
ENDIF()

IF (NEKTAR_BUILD_UNIT_TESTS)
    ADD_NEKPY_EXECUTABLE(test_nekmesh_node.py test_nekmesh_node.py)
    ADD_NEKPY_EXECUTABLE(test_nekmesh_element.py test_nekmesh_element.py)
    ADD_NEKPY_EXECUTABLE(test_nekmesh_mesh.py test_nekmesh_mesh.py)
    ADD_NEKPY_EXECUTABLE(test_nekmesh_module.py test_nekmesh_module.py)

    ADD_NEKPY_TEST(test_nekmesh_node)
    ADD_NEKPY_TEST(test_nekmesh_element)
    ADD_NEKPY_TEST(test_nekmesh_mesh)
    ADD_NEKPY_TEST(test_nekmesh_module)
ENDIF()
