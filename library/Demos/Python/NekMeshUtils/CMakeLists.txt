ADD_NEKPY_EXECUTABLE(test_nekmeshutil_node.py test_nekmeshutil_node.py)

IF (NEKTAR_BUILD_TESTS)
  ADD_NEKTAR_TEST(test_nekmeshutil_node)
ENDIF()
