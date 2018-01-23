SET(STDREGIONS_SOURCES
    StdRegions.cpp
    StdExpansion.cpp
    StdSegExp.cpp
    StdTriExp.cpp
    StdQuadExp.cpp
    )

ADD_NEKPY_LIBRARY(StdRegions STDREGIONS_SOURCES)
