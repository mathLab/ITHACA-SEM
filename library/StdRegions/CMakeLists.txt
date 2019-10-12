SET(STD_REGIONS_SOURCES
StdExpansion.cpp
StdExpansion0D.cpp
StdExpansion1D.cpp
StdExpansion2D.cpp
StdExpansion3D.cpp
StdHexExp.cpp
StdMatrixKey.cpp
StdNodalPrismExp.cpp
StdNodalTetExp.cpp
StdNodalTriExp.cpp
StdPrismExp.cpp
StdPyrExp.cpp
StdQuadExp.cpp
StdSegExp.cpp
StdPointExp.cpp
StdTetExp.cpp
StdTriExp.cpp
IndexMapKey.cpp
)

SET(STD_REGIONS_HEADERS
StdExpansion0D.h
StdExpansion1D.h
StdExpansion2D.h
StdExpansion3D.h
StdExpansion.h
StdHexExp.h
StdLinSysKey.hpp
StdMatrixKey.h
StdNodalPrismExp.h
StdNodalTetExp.h
StdNodalTriExp.h
StdPrismExp.h
StdPyrExp.h
StdQuadExp.h
StdRegionsDeclspec.h
StdSegExp.h
StdPointExp.h
StdTetExp.h
StdTriExp.h
StdRegions.hpp
IndexMapKey.h
)

ADD_DEFINITIONS(-DSTD_REGIONS_EXPORTS)

ADD_NEKTAR_LIBRARY(StdRegions
    SOURCES ${STD_REGIONS_SOURCES}
    HEADERS ${STD_REGIONS_HEADERS}
    DEPENDS LibUtilities
    SUMMARY "Nektar++ StdRegions library"
    DESCRIPTION "This library provides construction of the reference expansions for the various 1D, 2D and 3D regions.")

IF (NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/StdRegions COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

