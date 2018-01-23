SET(LOCALREGIONS_SOURCES
    LocalRegions.cpp
    Expansion.cpp
    QuadExp.cpp
    SegExp.cpp
    TriExp.cpp
    )

ADD_NEKPY_LIBRARY(LocalRegions LOCALREGIONS_SOURCES)
