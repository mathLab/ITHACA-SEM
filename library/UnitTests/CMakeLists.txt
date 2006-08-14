
SET(UnitTestSources
main.cpp
testNekManager.cpp
testNekMatrix.cpp
testNekVector.cpp
testNekPoint.cpp
)

SET(UnitTestHeaders
testNekManager.h
testNekMatrix.h
testNekPoint.h
testNekVector.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)

ADD_EXECUTABLE(UnitTests ${UnitTestSources} ${UnitTestHeaders})

TARGET_LINK_LIBRARIES(UnitTests
    optimized LibUtilities debug LibUtilities-g
    optimized ${BoostThreadName} debug ${BoostThreadDebugName}
)

INSTALL_TARGETS(${NEKTAR_BIN_DIR} UnitTests)