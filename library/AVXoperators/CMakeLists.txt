SET(AVXOPERATORS_SOURCES
    AVXIProduct.cpp
    AVXBwdTrans.cpp
    AVXPhysDeriv.cpp
    AVXHelmholtz.cpp
    AVXHelmholtzGlobal.cpp
    Operator.cpp
  )

SET(AVXOPERATORS_HEADERS
    AVXIProduct.h
    AVXIProductKernels.hpp
    AVXBwdTrans.h
    AVXBwdTransKernel.hpp
    AVXPhysDeriv.h
    AVXPhysDerivKernel.hpp
    AVXHelmholtz.h
    AVXHelmholtzGlobal.h
    AVXUtil.hpp
    VecData.hpp
    Operator.hpp
)

ADD_NEKTAR_LIBRARY(AVXoperators
    SOURCES ${AVXOPERATORS_SOURCES}
    HEADERS ${AVXOPERATORS_HEADERS}
    DEPENDS LibUtilities
    SUMMARY "Nektar++ AVX operators library"
    DESCRIPTION "This library provides key operators with AVX vector instructions.")

INSTALL(DIRECTORY ./
    DESTINATION ${NEKTAR_INCLUDE_DIR}/AVXoperators
    COMPONENT dev
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
