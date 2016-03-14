########################################################################
#
# ThirdParty configuration for Nektar++
#
# ANN
#
########################################################################



IF(NEKTAR_USE_MESH)

find_library(LIB_ANN ANN
    PATHS /opt/local/lib)

FIND_PATH(ANN_INCLUDE_DIR ANN.h
    HINTS /opt/local/include/ANN/
)

FIND_PATH(ANN_LIB_DIR libANN.a
    HINTS /opt/local/lib/
)

LINK_DIRECTORIES(${ANN_LIB_DIR})

INCLUDE_DIRECTORIES(SYSTEM ${ANN_INCLUDE_DIR})

ENDIF(NEKTAR_USE_MESH)
