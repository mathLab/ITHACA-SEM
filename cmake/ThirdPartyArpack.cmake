# Arpack
SET(NEKTAR_USE_ARPACK OFF CACHE BOOL 
    "Use Arpack routines for evaluating the eigenvalues and eigenvectors")

IF( NEKTAR_USE_ARPACK )
        INCLUDE (FindArpack)
        INCLUDE_DIRECTORIES(${ARPACK_INCLUDE_DIR})
        ADD_DEFINITIONS(-DNEKTAR_USING_ARPACK)
ENDIF( NEKTAR_USE_ARPACK)

