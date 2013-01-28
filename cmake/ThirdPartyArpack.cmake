# Arpack
OPTION(NEKTAR_USE_ARPACK
    "Use Arpack routines for evaluating the eigenvalues and eigenvectors" OFF)

IF( NEKTAR_USE_ARPACK )
    INCLUDE (FindArpack)
    INCLUDE_DIRECTORIES(${ARPACK_INCLUDE_DIR})
    ADD_DEFINITIONS(-DNEKTAR_USING_ARPACK)
ENDIF( NEKTAR_USE_ARPACK)

