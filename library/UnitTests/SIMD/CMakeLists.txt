SET(VecDataTestSources
    main.cpp
    TestVecData.cpp
    TestAVXUtil.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(VecDataUnitTests
    COMPONENT unit-test DEPENDS AVXOperators SOURCES ${VecDataTestSources})
ADD_TEST(NAME VecDataUnitTests COMMAND VecDataUnitTests --detect_memory_leaks=0)
