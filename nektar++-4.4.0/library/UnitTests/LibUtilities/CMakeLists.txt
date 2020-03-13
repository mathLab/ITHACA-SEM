
SET(PrecompiledHeaderSources
    TestConsistentObjectAccess.cpp
    TestLowerTriangularMatrix.cpp
    TestMatrixStoragePolicies.cpp
    TestNekMatrixMultiplication.cpp
    TestNekMatrixOperations.cpp
    TestRawType.cpp
    TestUpperTriangularMatrix.cpp
    TestSharedArray.cpp
    ../util.cpp
)

SET(UnitTestSources ${PrecompiledHeaderSources} main.cpp)   

SET(UnitTestHeaders
	LibUtilitiesUnitTestsPrecompiledHeader.h
	../util.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

ADD_NEKTAR_EXECUTABLE(LibUtilitiesUnitTests unit-test UnitTestSources UnitTestHeaders)
TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests LibUtilities)

#SET(PrecompiledHeaderName LibUtilitiesUnitTestsPrecompiledHeader.h)
#SETUP_PRECOMPILED_HEADERS(PrecompiledHeaderSources PrecompiledHeaderName)

ADD_TEST(NAME LibUtilities COMMAND LibUtilitiesUnitTests --detect_memory_leaks=0)
SUBDIRS(LinearAlgebra)

IF( NEKTAR_USE_EXPRESSION_TEMPLATES )
    SUBDIRS(ExpressionTemplates)
ENDIF()

