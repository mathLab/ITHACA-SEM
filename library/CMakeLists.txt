SET(LibrarySubDirs LibUtilities LocalRegions Collections
    MultiRegions SpatialDomains StdRegions SolverUtils)
SET(UnitTestSubDirs UnitTests)
SET(DemoSubDirs     Demos)
SET(TimingsSubDirs  Timings)

SUBDIRS(${LibrarySubDirs} )

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

IF( NEKTAR_BUILD_UNIT_TESTS )
    SUBDIRS(${UnitTestSubDirs})
ENDIF( NEKTAR_BUILD_UNIT_TESTS )

IF( NEKTAR_BUILD_DEMOS )
    SUBDIRS(${DemoSubDirs})
ENDIF( NEKTAR_BUILD_DEMOS )

IF( NEKTAR_BUILD_TIMINGS ) 
    SUBDIRS(${TimingsSubDirs})
ENDIF( NEKTAR_BUILD_TIMINGS ) 

