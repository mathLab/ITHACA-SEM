SET(LocalRegionsUnitTestsSources
    main.cpp
    TestGetCoords.cpp
    TestReOrientFacePhysMap.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(LocalRegionsUnitTests
    COMPONENT unit-test DEPENDS LocalRegions SOURCES ${LocalRegionsUnitTestsSources})
ADD_TEST(NAME LocalRegionsUnitTests COMMAND LocalRegionsUnitTests --detect_memory_leaks=0)
