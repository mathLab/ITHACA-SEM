#
# NektarSIMD.cmake
#
# Sets up cmake variables needed for the SIMD library in Nektar++
#

OPTION(NEKTAR_ENABLE_SIMD_SSE2 "Enable sse2 vector types" OFF)
OPTION(NEKTAR_ENABLE_SIMD_AVX2 "Enable avx2 vector types" OFF)
OPTION(NEKTAR_ENABLE_SIMD_AVX512 "Enable avx512 vector types" OFF)
MARK_AS_ADVANCED(FORCE NEKTAR_ENABLE_SIMD_SSE2 NEKTAR_ENABLE_SIMD_AVX2 NEKTAR_ENABLE_SIMD_AVX512)
IF (NEKTAR_ENABLE_SIMD_AVX512)
    MESSAGE(STATUS "Enabling avx512, you might need to clear CMAKE_CXX_FLAGS or add the appriopriate flags")
    ADD_DEFINITIONS(-DNEKTAR_ENABLE_SIMD_AVX512)
    SET(CMAKE_CXX_FLAGS "-mavx512 -mfma" CACHE STRING
        "Flags used by the CXX compiler during all build types.")
    SET(NEKTAR_ENABLE_SIMD_AVX2 "ON" FORCE)
    SET(NEKTAR_ENABLE_SIMD_SSE2 "ON" FORCE)
ENDIF()
IF (NEKTAR_ENABLE_SIMD_AVX2)
    MESSAGE(STATUS "Enabling avx2, you might need to clear CMAKE_CXX_FLAGS or add the appriopriate flags")
    ADD_DEFINITIONS(-DNEKTAR_ENABLE_SIMD_AVX2)
    SET(CMAKE_CXX_FLAGS "-mavx2 -mfma" CACHE STRING
        "Flags used by the CXX compiler during all build types.")
    SET(NEKTAR_ENABLE_SIMD_SSE2 "ON" FORCE)
ENDIF()
IF (NEKTAR_ENABLE_SIMD_SSE2)
    MESSAGE(STATUS "Enabling sse2, you might need to clear CMAKE_CXX_FLAGS or add the appriopriate flags")
    SET(SSE_FLAGS "-msse2")
    SET(CMAKE_CXX_FLAGS "${SSE_FLAGS}" CACHE STRING
        "Flags used by the CXX compiler during all build types.")
    ADD_DEFINITIONS(-DNEKTAR_ENABLE_SIMD_SSE2)
ENDIF()
# Vmath Simd
OPTION(NEKTAR_ENABLE_SIMD_VMATH "Enable vector types in vmath" OFF)
IF (NEKTAR_ENABLE_SIMD_VMATH)
    ADD_DEFINITIONS(-DNEKTAR_ENABLE_SIMD_VMATH)
ENDIF()
MARK_AS_ADVANCED(FORCE NEKTAR_ENABLE_SIMD_VMATH)