SET(LibUtilitiesUnitTestSources
    TestConsistentObjectAccess.cpp
    TestRawType.cpp
    TestInterpreter.cpp
    TestSharedArray.cpp
    TestCheckedCast.cpp
    TestRealComparison.cpp
    main.cpp
    ../util.cpp
)

SET(UnitTestSources ${PrecompiledHeaderSources} main.cpp)

SET(UnitTestHeaders
	LibUtilitiesUnitTestsPrecompiledHeader.h
	../util.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(LibUtilitiesUnitTests
    COMPONENT unit-test DEPENDS LibUtilities SOURCES ${LibUtilitiesUnitTestSources})
ADD_TEST(NAME LibUtilitiesUnitTests COMMAND LibUtilitiesUnitTests --detect_memory_leaks=0)

SUBDIRS(LinearAlgebra)
