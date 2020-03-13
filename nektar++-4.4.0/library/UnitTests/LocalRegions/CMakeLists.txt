
SET(Sources
    main.cpp
    TestGetCoords.cpp
)

SET(Headers
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

SET(ProjectName LocalRegionsUnitTests)
ADD_NEKTAR_EXECUTABLE(${ProjectName} unit-test Sources Headers)

TARGET_LINK_LIBRARIES(${ProjectName} LocalRegions)

ADD_TEST(NAME LocalRegions COMMAND ${ProjectName} --detect_memory_leaks=0)
