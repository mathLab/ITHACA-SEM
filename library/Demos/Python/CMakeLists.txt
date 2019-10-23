ADD_NEKPY_EXECUTABLE(StdIntegral2D.py StdIntegral2D.py)
ADD_NEKPY_EXECUTABLE(Helmholtz2D.py Helmholtz2D.py)
ADD_NEKPY_EXECUTABLE(CurvedGeom.py CurvedGeom.py)
ADD_NEKPY_EXECUTABLE(NekPy_ReferenceTest_CppDeleteFirst.py NekPy_ReferenceTest_CppDeleteFirst.py)
ADD_NEKPY_EXECUTABLE(NekPy_ReferenceTest_PythonDeleteFirst.py NekPy_ReferenceTest_PythonDeleteFirst.py)
ADD_NEKPY_EXECUTABLE(NekPy_MemoryTest_PythonToCpp_Address.py NekPy_MemoryTest_PythonToCpp_Address.py)

IF (NEKTAR_BUILD_TESTS)
    ADD_NEKTAR_TEST(StdIntegral2D)

    ADD_NEKTAR_TEST(CurvedGeom)
    ADD_NEKTAR_TEST(NekPy_ReferenceTest_CppDeleteFirst)
    ADD_NEKTAR_TEST(NekPy_ReferenceTest_PythonDeleteFirst)
    ADD_NEKTAR_TEST(NekPy_MemoryTest_PythonToCpp_Address)

    IF (NEKTAR_USE_MPI)
        ADD_NEKTAR_TEST(Helmholtz2D)
    ENDIF()
ENDIF()
