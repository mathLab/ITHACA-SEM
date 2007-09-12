
SET(UnitTestSources
    CountedObject.cpp
    main.cpp
    ScaledMatrixUnitTests.cpp
    BlockMatrixUnitTests.cpp
    testBoostUtil.cpp
    testExpressionTemplates.cpp
    testLinearSystem.cpp
    testNekLinAlgAlgorithms.cpp
    testNekManager.cpp
    testNekMatrix.cpp
    testNekMatrixOperations.cpp
    testNekPoint.cpp
    TestNekPtr.cpp
    testNekVector.cpp
    testNekSharedArray.cpp
    Memory/TestNekMemoryManager.cpp
    StdRegions/testStdSegExp.cpp
    testFoundation/testFoundation.cpp
    testFoundation/testInterpolation.cpp
    testFoundation/testDerivation.cpp
)

SET(UnitTestHeaders
    CountedObject.h
    testBoostUtil.h
    testExpressionTemplates.h
    testLinearSystem.h
    testNekLinAlgAlgorithms.h
    testNekManager.h
    testNekMatrix.h
    testNekMatrixOperations.h
    testNekPoint.h
    TestNekPtr.h
    testNekVector.h
    testNekSharedArray.h
    Memory/TestNekMemoryManager.h
    StdRegions/testStdSegExp.h
    testFoundation/testFoundation.h
    testFoundation/testInterpolation.h
    testFoundation/testDerivation.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)

ADD_NEKTAR_EXECUTABLE(UnitTests UnitTestSources UnitTestHeaders)

IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF( USE_GCOV )
		TARGET_LINK_LIBRARIES(UnitTests gcov)
	ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

TARGET_LINK_LIBRARIES(UnitTests
    optimized LibUtilities debug LibUtilities-g
    optimized StdRegions debug StdRegions-g
    optimized ${BOOST_THREAD_LIB} debug ${BOOST_THREAD_DEBUG_LIB}
)

SET_LAPACK_LINK_LIBRARIES(UnitTests)

IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF( USE_GCOV )
        TARGET_LINK_LIBRARIES(UnitTests gcov)
    ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

SUBDIRS(LibUtilities)
INSTALL(TARGETS UnitTests RUNTIME DESTINATION ${NEKTAR_BIN_DIR})
