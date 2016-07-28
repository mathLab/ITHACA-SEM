SET(LibrarySubDirs FieldUtils GlobalMapping LibUtilities LocalRegions
    Collections MultiRegions SpatialDomains StdRegions)
SET(UnitTestSubDirs UnitTests)
SET(DemoSubDirs     Demos)
SET(TimingsSubDirs  Timings)
SET(SolverUtilsSubDirs SolverUtils)
SET(NekMeshUtilsSubDirs NekMeshUtils)

SUBDIRS(${LibrarySubDirs})

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

IF( NEKTAR_BUILD_SOLVERS )
    SUBDIRS(${SolverUtilsSubDirs})
ENDIF( NEKTAR_BUILD_SOLVERS )

IF( NEKTAR_BUILD_UTILITIES )
    SUBDIRS(${NekMeshUtilsSubDirs})
ENDIF( NEKTAR_BUILD_UTILITIES )
