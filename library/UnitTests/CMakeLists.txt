
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

SET(UnitTestHeaders
    CountedObject.h
    testLinearSystem.h
    testNekLinAlgAlgorithms.h
    testNekManager.h
    testNekPoint.h
    testNekSharedArray.h
    Memory/TestNekMemoryManager.h
    StdRegions/testStdSegExp.h
    testFoundation/testFoundation.h
    testFoundation/testInterpolation.h
    testFoundation/testDerivation.h
    util.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS -DNEKTAR_UNIT_TESTS)

ADD_NEKTAR_EXECUTABLE(UnitTests unit-test UnitTestSources UnitTestHeaders)

TARGET_LINK_LIBRARIES(UnitTests LibUtilities StdRegions)

#ADD_TEST(NAME UnitTests COMMAND UnitTests --detect_memory_leaks=0)
SUBDIRS(LibUtilities LocalRegions)
