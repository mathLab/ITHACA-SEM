# Try to find OCE

set(TEST_ENV $ENV{OCE_ROOT})
set(TEST_ENV1 $ENV{OCE_DIR})
if(NOT DEFINED TEST_ENV AND NOT DEFINED TEST_ENV1)
    message(STATUS "manually scan for oce")
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(ENV{OCE_DIR} "/opt/local/Library/Frameworks")
    elseif(UNIX)
        set(ENV{OCE_DIR} "/usr/local/share/cmake/")
    endif()
endif()

find_package(OCE 0.17 QUIET)

if(OCE_FOUND)
  message(STATUS "Found OpenCASCADE Community Edition. Version ${OCE_VERSION}")

  file(STRINGS ${OCE_INCLUDE_DIRS}/Standard_Version.hxx OCC_MAJOR
    REGEX "#define OCC_VERSION_MAJOR.*"
  )
  string(REGEX MATCH "[0-9]+" OCC_MAJOR ${OCC_MAJOR})

  file(STRINGS ${OCE_INCLUDE_DIRS}/Standard_Version.hxx OCC_MINOR
    REGEX "#define OCC_VERSION_MINOR.*"
  )
  string(REGEX MATCH "[0-9]+" OCC_MINOR ${OCC_MINOR})

  file(STRINGS ${OCE_INCLUDE_DIRS}/Standard_Version.hxx OCC_MAINT
    REGEX "#define OCC_VERSION_MAINTENANCE.*"
  )
  string(REGEX MATCH "[0-9]+" OCC_MAINT ${OCC_MAINT})

  set(OCC_VERSION_STRING "${OCC_MAJOR}.${OCC_MINOR}.${OCC_MAINT}")

  message(STATUS "-- OCC version: ${OCC_VERSION_STRING}")

  set(OCE_LIBRARIES
    TKFillet
    TKMesh
    TKernel
    TKG2d
    TKG3d
    TKMath
    TKIGES
    TKSTL
    TKShHealing
    TKXSBase
    TKBool
    TKBO
    TKBRep
    TKTopAlgo
    TKGeomAlgo
    TKGeomBase
    TKOffset
    TKPrim
    TKSTEP
    TKSTEPBase
    TKSTEPAttr
    TKHLR
    TKFeat
  )
endif(OCE_FOUND)
