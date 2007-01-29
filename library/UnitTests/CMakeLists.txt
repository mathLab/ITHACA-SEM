
SET(UnitTestSources
CountedObject.cpp
main.cpp
testBoostUtil.cpp
testExpressionTemplates.cpp
testLinearSystem.cpp
testNekLinAlgAlgorithms.cpp
testNekManager.cpp
testNekMatrix.cpp
testNekPoint.cpp
testNekVector.cpp
Memory/TestNekMemoryManager.cpp
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
Memory/TestNekMemoryManager.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)

ADD_EXECUTABLE(UnitTests ${UnitTestSources} ${UnitTestHeaders})

IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF( USE_GCOV )
		TARGET_LINK_LIBRARIES(UnitTests gcov)
	ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

TARGET_LINK_LIBRARIES(UnitTests
    optimized LibUtilities debug LibUtilities-g
    optimized ${BoostThreadName} debug ${BoostThreadDebugName}
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

INSTALL_TARGETS(${NEKTAR_BIN_DIR} UnitTests)
