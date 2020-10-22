ADD_NEKPY_EXECUTABLE(ReferenceTest_CppDeleteFirst.py ReferenceTest_CppDeleteFirst.py)
ADD_NEKPY_EXECUTABLE(ReferenceTest_PythonDeleteFirst.py ReferenceTest_PythonDeleteFirst.py)
ADD_NEKPY_EXECUTABLE(MemoryTest_PythonToCpp_Address.py MemoryTest_PythonToCpp_Address.py)

IF (NEKTAR_BUILD_TESTS)
    ADD_NEKPY_TEST(ReferenceTest_CppDeleteFirst)
    ADD_NEKPY_TEST(ReferenceTest_PythonDeleteFirst)
    ADD_NEKPY_TEST(MemoryTest_PythonToCpp_Address)
ENDIF()
