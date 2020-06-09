ADD_NEKPY_EXECUTABLE(test_nekmesh_node.py test_nekmesh_node.py)
ADD_NEKPY_EXECUTABLE(test_nekmesh_element.py test_nekmesh_element.py)
ADD_NEKPY_EXECUTABLE(test_nekmesh_mesh.py test_nekmesh_mesh.py)

IF (NEKTAR_BUILD_TESTS)
  ADD_NEKTAR_TEST(test_nekmesh_node)
  ADD_NEKTAR_TEST(test_nekmesh_element)
  ADD_NEKTAR_TEST(test_nekmesh_mesh)
ENDIF()
