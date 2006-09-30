
SET(UnitTestSources
main.cpp
testNekManager.cpp
testNekMatrix.cpp
testNekVector.cpp
testNekPoint.cpp
testExpressionTemplates.cpp
testLinearSystem.cpp
)

SET(UnitTestHeaders
testNekManager.h
testNekMatrix.h
testNekPoint.h
testNekVector.h
testExpressionTemplates.h
testLinearSystem.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS -DNEKTAR_USE_EXPRESSION_TEMPLATES)

ADD_EXECUTABLE(UnitTests ${UnitTestSources} ${UnitTestHeaders})

TARGET_LINK_LIBRARIES(UnitTests
    optimized LibUtilities debug LibUtilities-g
    optimized ${BoostThreadName} debug ${BoostThreadDebugName}
)

IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF( USE_GCOV )
        TARGET_LINK_LIBRARIES(UnitTests gcov)
    ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

INSTALL_TARGETS(${NEKTAR_BIN_DIR} UnitTests)
