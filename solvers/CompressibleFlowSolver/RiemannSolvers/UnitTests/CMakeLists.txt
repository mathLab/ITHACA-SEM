SET(RiemannCompUnitTestsSources
    main.cpp
    TestRiemann.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(RiemannCompUnitTests
    COMPONENT unit-test
    DEPENDS CompRiemann
    SOURCES ${RiemannCompUnitTestsSources})
ADD_TEST(NAME RiemannCompUnitTests COMMAND RiemannCompUnitTests --detect_memory_leaks=0)