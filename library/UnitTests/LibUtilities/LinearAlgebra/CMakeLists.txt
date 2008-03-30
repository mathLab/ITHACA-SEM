
SET(TARGET_NAME LinearAlgebraUnitTests)

SET(Sources
    main.cpp   
    TestBandedMatrixOperations.cpp
    TestBandedMatrixStoragePolicy.cpp
    TestBlockMatrix.cpp
    TestDiagonalMatrixStoragePolicy.cpp
    TestDiagonalMatrixOperations.cpp
    TestFullMatrixStoragePolicy.cpp
    TestFullMatrixOperations.cpp
    TestLowerTriangularMatrixStoragePolicy.cpp
    TestNekVector.cpp
    TestScaledBlockMatrixOperations.cpp
    TestScaledMatrix.cpp
    TestSymmetricMatrixStoragePolicy.cpp
    TestUpperTriangularMatrixStoragePolicy.cpp
    TestStandardMatrix.cpp
)

SET(Headers
    TestCombinationRunner.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)

ADD_NEKTAR_EXECUTABLE(${TARGET_NAME} Sources Headers)


TARGET_LINK_LIBRARIES(${TARGET_NAME}
    optimized LibUtilities debug LibUtilities-g
    optimized ${BOOST_THREAD_LIB} debug ${BOOST_THREAD_DEBUG_LIB}
)

SET_LAPACK_LINK_LIBRARIES(${TARGET_NAME})

