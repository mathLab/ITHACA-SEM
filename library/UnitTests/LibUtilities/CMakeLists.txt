SET(LibUtilitiesUnitTestSources
    TestConsistentObjectAccess.cpp
    TestLowerTriangularMatrix.cpp
    TestMatrixStoragePolicies.cpp
    TestNekMatrixMultiplication.cpp
    TestNekMatrixOperations.cpp
    TestRawType.cpp
    TestUpperTriangularMatrix.cpp
    TestSharedArray.cpp
    main.cpp
    ../util.cpp
)

SET(UnitTestSources ${PrecompiledHeaderSources} main.cpp)   

SET(UnitTestHeaders
	LibUtilitiesUnitTestsPrecompiledHeader.h
	../util.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(LibUtilitiesUnitTests
    COMPONENT unit-test DEPENDS LibUtilities SOURCES ${LibUtilitiesUnitTestSources})
ADD_TEST(NAME LibUtilitiesUnitTests COMMAND LibUtilitiesUnitTests --detect_memory_leaks=0)

SUBDIRS(LinearAlgebra)

IF( NEKTAR_USE_EXPRESSION_TEMPLATES )
    SUBDIRS(ExpressionTemplates)
ENDIF()

