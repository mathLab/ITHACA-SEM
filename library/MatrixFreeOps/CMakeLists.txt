SET(MATRIXFREEOPS_SOURCES
    BwdTrans.cpp
    Helmholtz.cpp
    IProduct.cpp
    IProductWRTDerivBase.cpp
    PhysDeriv.cpp
    # AVXHelmholtzGlobal.cpp
    Operator.cpp
  )

SET(MATRIXFREEOPS_HEADERS
    # AVXAssembly.h
    BwdTrans.h
    BwdTransKernels.hpp
    Helmholtz.h
    IProduct.h
    IProductKernels.hpp
    IProductWRTDerivBase.h
    # AVXHelmholtzGlobal.h
    Util.hpp
    Operator.hpp
)

ADD_NEKTAR_LIBRARY(MatrixFreeOps
    SOURCES ${MATRIXFREEOPS_SOURCES}
    HEADERS ${MATRIXFREEOPS_HEADERS}
    DEPENDS LibUtilities
    SUMMARY "Nektar++ matrix free operators library"
    DESCRIPTION "This library provides key matrix free operators that can use vector instructions.")


INSTALL(DIRECTORY ./
    DESTINATION ${NEKTAR_INCLUDE_DIR}/MatrixFreeOps
    COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
