
SET(Sources
    main.cpp
    TestGetCoords.cpp
)

SET(Headers
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

SET(ProjectName LocalRegionsUnitTests)
ADD_NEKTAR_EXECUTABLE(${ProjectName} Sources Headers)

TARGET_LINK_LIBRARIES(${ProjectName}
    LibUtilities
    StdRegions
    LocalRegions
    SpatialDomains
    ${Boost_THREAD_LIBRARY}
)

SET_LAPACK_LINK_LIBRARIES(${ProjectName})

