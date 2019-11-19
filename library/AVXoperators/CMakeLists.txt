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
    AVXBwdTransKernels.hpp
    AVXPhysDeriv.h
    AVXPhysDerivKernels.hpp
    AVXHelmholtz.h
    AVXHelmholtzGlobal.h
    AVXUtil.hpp
    Operator.hpp
)

ADD_NEKTAR_LIBRARY(AVXoperators
    SOURCES ${AVXOPERATORS_SOURCES}
    HEADERS ${AVXOPERATORS_HEADERS}
    DEPENDS LocalRegions
    SUMMARY "Nektar++ AVX operators library"
    DESCRIPTION "This library provides key operators with AVX vector instructions.")


TARGET_LINK_LIBRARIES(AVXoperators
    ${NEKTAR++_LIBRARIES}
    ${NEKTAR++_TP_LIBRARIES}
    ${LIBXSMM_LIBRARY}
    ${ASMJIT_LIBRARY})

IF (LIKWID_LIBRARY)
    TARGET_LINK_LIBRARIES(AVXoperators ${LIKWID_LIBRARY})
ENDIF()

INSTALL(DIRECTORY ./
    DESTINATION ${NEKTAR_INCLUDE_DIR}/AVXoperators
    COMPONENT dev
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
