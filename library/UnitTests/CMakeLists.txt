SET(UnitTestSources
    main.cpp
    ScaledMatrixUnitTests.cpp
    BlockMatrixUnitTests.cpp
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
    util.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS -DNEKTAR_UNIT_TESTS)

ADD_NEKTAR_EXECUTABLE(UnitTests
    COMPONENT unit-test DEPENDS StdRegions SOURCES ${UnitTestSources})

ADD_TEST(NAME UnitTests COMMAND UnitTests --detect_memory_leaks=0)
SUBDIRS(LibUtilities LocalRegions Collections)