
SET(TARGET_NAME LinearAlgebraUnitTests)

SET(Sources
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
    ../../util.cpp
)

SET(Headers
    TestCombinationRunner.h
    ../../util.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)

ADD_NEKTAR_EXECUTABLE(${TARGET_NAME} unit-test Sources Headers)
TARGET_LINK_LIBRARIES(${TARGET_NAME} LibUtilities)

#ADD_TEST(NAME LinearAlgebra COMMAND ${TARGET_NAME} --detect_memory_leaks=0)
