
SET(UnitTestSources
main.cpp
testNekMatrix.cpp
testNekVector.cpp
testNekPoint.cpp
)

SET(UnitTestHeaders
testNekMatrix.h
testNekPoint.h
testNekVector.h
)

ADD_EXECUTABLE(UnitTests ${UnitTestSources} ${UnitTestsHeaders})

TARGET_LINK_LIBRARIES(UnitTests
    optimized LibUtilities debug LibUtilities-g
    optimized boost_thread-gcc-mt debug boost_thread-gcc-mt-d
)

INSTALL_TARGETS(${NEKTAR_BIN_DIR} UnitTests)