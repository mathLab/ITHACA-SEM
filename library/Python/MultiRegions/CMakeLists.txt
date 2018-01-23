SET(MULTIREGIONS_SOURCES
    MultiRegions.cpp
    ExpList.cpp
    ExpList2D.cpp
    )

ADD_NEKPY_LIBRARY(MultiRegions MULTIREGIONS_SOURCES)
