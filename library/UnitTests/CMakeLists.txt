
SET(UnitTestSources
CountedObject.cpp
main.cpp
ScaledMatrixUnitTests.cpp
testBoostUtil.cpp
testExpressionTemplates.cpp
testLinearSystem.cpp
testNekLinAlgAlgorithms.cpp
testNekManager.cpp
testNekMatrix.cpp
testNekPoint.cpp
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
testNekPoint.h
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

#ADD_NEKTAR_EXECUTABLE(UnitTests UnitTestSources)

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


IF( USE_MKL AND MKL_FOUND )
        TARGET_LINK_LIBRARIES(UnitTests
            ${MKL_LAPACK}
            optimized ${MKL} debug ${MKL}
            ${MKL_GUIDE}
        )
ENDIF( USE_MKL AND MKL_FOUND )

IF( USE_ACML AND ACML_FOUND )
    TARGET_LINK_LIBRARIES(UnitTests ${ACML_TARGET_LINK_LIBRARIES}  )
ENDIF( USE_ACML AND ACML_FOUND )

IF( USE_ATLAS AND ATLAS_FOUND )
        TARGET_LINK_LIBRARIES(UnitTests
            ${ATLAS_LAPACK}
            ${ATLAS_CBLAS}
            ${ATLAS} )
ENDIF( USE_ATLAS AND ATLAS_FOUND )

IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF( USE_GCOV )
        TARGET_LINK_LIBRARIES(UnitTests gcov)
    ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

IF( USE_FRAMEWORK )
    TARGET_LINK_LIBRARIES(UnitTests ${FRAMEWORK_LINK_FLAGS})
ENDIF ( USE_FRAMEWORK )

IF( USE_WIN32_LAPACK )
	TARGET_LINK_LIBRARIES(UnitTests ${WIN32_LAPACK} ${WIN32_BLAS})	
ENDIF( USE_WIN32_LAPACK )

INSTALL(TARGETS UnitTests RUNTIME DESTINATION ${NEKTAR_BIN_DIR})
