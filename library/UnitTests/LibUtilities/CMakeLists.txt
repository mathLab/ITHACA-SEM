
SET(PrecompiledHeaderSources
	TestLowerTriangularMatrix.cpp
    TestMatrixStoragePolicies.cpp
    TestNekMatrixMultiplication.cpp
    TestNekMatrixOperations.cpp
    TestUpperTriangularMatrix.cpp
)

SET(UnitTestSources ${PrecompiledHeaderSources} main.cpp)   

SET(UnitTestHeaders
	LibUtilitiesUnitTestsPrecompiledHeader.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(LibUtilitiesUnitTests UnitTestSources UnitTestHeaders)

SET(PrecompiledHeaderName LibUtilitiesUnitTestsPrecompiledHeader.h)
SETUP_PRECOMPILED_HEADERS(PrecompiledHeaderSources PrecompiledHeaderName)


IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF( USE_GCOV )
		TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests gcov)
	ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests
    optimized LibUtilities debug LibUtilities-g
    optimized StdRegions debug StdRegions-g
    optimized ${BOOST_THREAD_LIB} debug ${BOOST_THREAD_DEBUG_LIB}
)

SET_LAPACK_LINK_LIBRARIES(LibUtilitiesUnitTests)

INSTALL(TARGETS LibUtilitiesUnitTests RUNTIME DESTINATION ${NEKTAR_BIN_DIR})

SUBDIRS(LinearAlgebra ExpressionTemplates)
