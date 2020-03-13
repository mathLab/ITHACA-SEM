########################################################################
#
# ThirdParty configuration for Nektar++
#
# ARPACK
#
########################################################################

OPTION(NEKTAR_USE_ARPACK
    "Use Arpack routines for evaluating eigenvalues and eigenvectors" OFF)

IF (NEKTAR_USE_ARPACK)
    FIND_LIBRARY(ARPACK_LIBRARY NAMES "arpack.1" "arpack" PATHS /opt/local/lib)

    IF (ARPACK_LIBRARY)
        MESSAGE(STATUS "Found Arpack: ${ARPACK_LIBRARY}")
        MARK_AS_ADVANCED(ARPACK_LIBRARY)
    ELSE()
        MESSAGE(FATAL_ERROR "Could not find Arpack")
    ENDIF()
ENDIF()

