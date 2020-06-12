SET(MATRIXFREEOPS_SOURCES
    AVXBwdTrans.cpp
    AVXIProduct.cpp
    AVXIProductWRTDerivBase.cpp
    AVXPhysDeriv.cpp
    # AVXHelmholtz.cpp
    # AVXHelmholtzGlobal.cpp
    Operator.cpp
  )

SET(MATRIXFREEOPS_HEADERS
    # AVXAssembly.h
    AVXBwdTrans.h
    AVXBwdTransKernels.hpp
    AVXIProduct.h
    AVXIProductKernels.hpp
    AVXIProductWRTDerivBase.h
    # AVXHelmholtzGlobal.h
    AVXUtil.hpp
    Operator.hpp
    # VecData.hpp
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
