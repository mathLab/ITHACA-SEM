SET(UnitTestSources
    main.cpp
    ScaledMatrixUnitTests.cpp
    BlockMatrixUnitTests.cpp
    testLinearSystem.cpp
    testNekLinAlgAlgorithms.cpp
    testNekMatrix.cpp
    testNekPoint.cpp
    testNekVector.cpp
    testNekSharedArray.cpp
    Memory/TestNekMemoryManager.cpp
    util.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS -DNEKTAR_UNIT_TESTS)

ADD_NEKTAR_EXECUTABLE(UnitTests
    COMPONENT unit-test DEPENDS StdRegions SOURCES ${UnitTestSources})

ADD_TEST(NAME UnitTests COMMAND UnitTests --detect_memory_leaks=0)
SUBDIRS(LibUtilities LocalRegions Collections)
