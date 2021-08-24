SUBDIRS(LibUtilities StdRegions SpatialDomains LocalRegions Collections MultiRegions)

IF(NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()
