SET(LinearAlgebraUnitTestsSources
    main.cpp   
    TestBandedMatrixOperations.cpp
    TestBandedMatrixStoragePolicy.cpp
    TestBlockMatrix.cpp
    TestCanGetRawPtr.cpp
    TestDgemmOptimizations.cpp
    TestDiagonalMatrixStoragePolicy.cpp
    TestDiagonalMatrixOperations.cpp
    TestFullMatrixStoragePolicy.cpp
    TestFullMatrixOperations.cpp
    TestLowerTriangularMatrixStoragePolicy.cpp
    TestNekVector.cpp
    TestScaledBlockMatrixOperations.cpp
    TestScaledMatrix.cpp
    TestSymmetricMatrixStoragePolicy.cpp
    TestTriangularMatrixOperations.cpp
    TestUpperTriangularMatrixStoragePolicy.cpp
    TestStandardMatrix.cpp
    TestLowerTriangularMatrix.cpp
    TestNekMatrixMultiplication.cpp
    TestNekMatrixOperations.cpp
    TestUpperTriangularMatrix.cpp
    ../../util.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(LinearAlgebraUnitTests
    COMPONENT unit-test DEPENDS LibUtilities SOURCES ${LinearAlgebraUnitTestsSources})
ADD_TEST(NAME LinearAlgebraUnitTests COMMAND LinearAlgebraUnitTests --detect_memory_leaks=0)
