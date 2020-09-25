ADD_NEKPY_EXECUTABLE(Helmholtz2D.py Helmholtz2D.py)
ADD_NEKPY_EXECUTABLE(CurvedGeom.py CurvedGeom.py)

IF (NEKTAR_BUILD_TESTS)
    ADD_NEKPY_TEST(CurvedGeom)

    IF (NEKTAR_USE_MPI)
        ADD_NEKPY_TEST(Helmholtz2D)
    ENDIF()
ENDIF()
