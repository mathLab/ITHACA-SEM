SET(CollectionsUnitTestsSources
    main.cpp
    # TestHexCollection.cpp
    TestQuadCollection.cpp
    # TestPrismCollection.cpp
    # TestPyrCollection.cpp
    TestSegCollection.cpp
    # TestTetCollection.cpp
    # TestTriCollection.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(CollectionsUnitTests
    COMPONENT unit-test DEPENDS Collections SOURCES ${CollectionsUnitTestsSources})
ADD_TEST(NAME CollectionsUnitTests COMMAND CollectionsUnitTests --detect_memory_leaks=0)
