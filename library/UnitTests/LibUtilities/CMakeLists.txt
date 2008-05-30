
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

ADD_NEKTAR_EXECUTABLE(LibUtilitiesUnitTests UnitTestSources UnitTestHeaders)

SET(PrecompiledHeaderName LibUtilitiesUnitTestsPrecompiledHeader.h)
SETUP_PRECOMPILED_HEADERS(PrecompiledHeaderSources PrecompiledHeaderName)

TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests
    optimized LibUtilities debug LibUtilities-g
    optimized StdRegions debug StdRegions-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
)

SET_LAPACK_LINK_LIBRARIES(LibUtilitiesUnitTests)


SUBDIRS(LinearAlgebra ExpressionTemplates)
