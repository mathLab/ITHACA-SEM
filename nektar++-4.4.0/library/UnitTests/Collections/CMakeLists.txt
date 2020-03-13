SET(Sources
    main.cpp
    TestHexCollection.cpp
    TestQuadCollection.cpp
    TestPrismCollection.cpp
    TestSegCollection.cpp
    TestTetCollection.cpp
    TestTriCollection.cpp
)

SET(Headers
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

SET(ProjectName CollectionsUnitTests)
ADD_NEKTAR_EXECUTABLE(${ProjectName} unit-test Sources Headers)

TARGET_LINK_LIBRARIES(${ProjectName} Collections)
                                     
ADD_TEST(NAME Collections COMMAND ${ProjectName} --detect_memory_leaks=0)
