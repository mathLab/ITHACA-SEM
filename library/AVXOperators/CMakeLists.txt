SET(AVXOPERATORS_SOURCES
    AVXBwdTrans.cpp
    AVXIProduct.cpp
    AVXIProductWRTDerivBase.cpp
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
    AVXIProductWRTDerivBase.h
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
