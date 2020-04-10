SET(AVXOPERATORS_SOURCES
    AVXBwdTrans.cpp
    AVXIProduct.cpp
    # AVXIProductDeriv.cpp
    AVXPhysDeriv.cpp
    # AVXHelmholtz.cpp
    # AVXHelmholtzGlobal.cpp
    Operator.cpp
  )

SET(AVXOPERATORS_HEADERS
    # AVXAssembly.h
    AVXBwdTrans.h
    AVXBwdTransKernels.hpp
    AVXIProduct.h
    AVXIProductKernels.hpp
    # AVXIProductDeriv.h
    # AVXIProductDerivKernels.hpp
    AVXPhysDeriv.h
    AVXPhysDerivKernels.hpp
    # AVXHelmholtz.h
    # AVXHelmholtzGlobal.h
    AVXUtil.hpp
    Operator.hpp
    VecData.hpp
)

ADD_NEKTAR_LIBRARY(AVXOperators
    SOURCES ${AVXOPERATORS_SOURCES}
    HEADERS ${AVXOPERATORS_HEADERS}
    DEPENDS LibUtilities
    SUMMARY "Nektar++ AVX operators library"
    DESCRIPTION "This library provides key operators with AVX vector instructions.")


INSTALL(DIRECTORY ./
    DESTINATION ${NEKTAR_INCLUDE_DIR}/AVXOperators
    COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
