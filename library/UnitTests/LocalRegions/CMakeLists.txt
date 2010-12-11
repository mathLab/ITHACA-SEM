
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
    optimized LibUtilities debug LibUtilities-g
    optimized StdRegions debug StdRegions-g
    optimized LocalRegions debug LocalRegions-g
    optimized SpatialDomains debug SpatialDomains-g
    ${Boost_THREAD_LIBRARY}
)

SET_LAPACK_LINK_LIBRARIES(${ProjectName})

