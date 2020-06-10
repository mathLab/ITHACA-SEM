SET(SIMDOPERATORS_SOURCES
    AVXBwdTrans.cpp
    AVXIProduct.cpp
    AVXIProductWRTDerivBase.cpp
    AVXPhysDeriv.cpp
    # AVXHelmholtz.cpp
    # AVXHelmholtzGlobal.cpp
    Operator.cpp
  )

SET(SIMDOPERATORS_HEADERS
    # AVXAssembly.h
    AVXBwdTrans.h
    AVXBwdTransKernels.hpp
    AVXIProduct.h
    AVXIProductKernels.hpp
    AVXIProductWRTDerivBase.h
    # AVXHelmholtzGlobal.h
    AVXUtil.hpp
    Operator.hpp
    VecData.hpp
)

ADD_NEKTAR_LIBRARY(SimdOperators
    SOURCES ${SIMDOPERATORS_SOURCES}
    HEADERS ${SIMDOPERATORS_HEADERS}
    DEPENDS LibUtilities
    SUMMARY "Nektar++ SIMD operators library"
    DESCRIPTION "This library provides key operators with SIMD vector instructions.")


INSTALL(DIRECTORY ./
    DESTINATION ${NEKTAR_INCLUDE_DIR}/SimdOperators
    COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
