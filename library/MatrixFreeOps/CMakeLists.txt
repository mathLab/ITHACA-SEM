SET(MATRIXFREEOPS_SOURCES
    BwdTrans.cpp
    Helmholtz.cpp
    IProductQuad.cpp
    IProductTri.cpp
    IProductTet.cpp
    IProductPrism.cpp
    IProductHex.cpp
    IProductWRTDerivBaseQuad.cpp
    IProductWRTDerivBaseTri.cpp
    IProductWRTDerivBaseTet.cpp
    IProductWRTDerivBasePrism.cpp
    IProductWRTDerivBaseHex.cpp
    PhysDeriv.cpp
    Operator.cpp
  )

SET(MATRIXFREEOPS_HEADERS
    BwdTrans.h
    BwdTransKernels.hpp
    Helmholtz.h
    IProduct.h
    IProductKernels.hpp
    IProductWRTDerivBase.h
    Operator.hpp
)

ADD_DEFINITIONS(-DMATRIXFREE_EXPORTS)

ADD_NEKTAR_LIBRARY(MatrixFreeOps
    SOURCES ${MATRIXFREEOPS_SOURCES}
    HEADERS ${MATRIXFREEOPS_HEADERS}
    DEPENDS LibUtilities
    SUMMARY "Nektar++ matrix free operators library"
    DESCRIPTION "This library provides key matrix free operators that can use vector instructions.")


INSTALL(DIRECTORY ./
    DESTINATION ${NEKTAR_INCLUDE_DIR}/MatrixFreeOps
    COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
