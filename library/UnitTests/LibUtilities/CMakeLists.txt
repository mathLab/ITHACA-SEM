
SET(UnitTestSources
    TestNekMatrixMultiplication.cpp
    TestNekMatrixOperations.cpp
	TestMatrixStoragePolicies.cpp
	TestUpperTriangularMatrix.cpp
    main.cpp   
)

SET(UnitTestHeaders
    TestNekMatrixMultiplication.h
    TestNekMatrixOperations.h
    TestMatrixStoragePolicies.h
    TestUpperTriangularMatrix.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)

ADD_NEKTAR_EXECUTABLE(LibUtilitiesUnitTests UnitTestSources UnitTestHeaders)


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
