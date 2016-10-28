# Try to find OCE

set(TEST_ENV $ENV{OCE_ROOT})
if(NOT DEFINED OCE_DIR AND DEFINED TEST_ENV)
  file(GLOB OCE_DIR $ENV{OCE_ROOT}/lib/oce-*)
endif()

# First try to find OpenCASCADE Community Edition
if(NOT DEFINED OCE_DIR)
  # Check for OSX needs to come first because UNIX evaluates to true on OSX
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if(DEFINED MACPORTS_PREFIX)
      find_package(OCE 0.17 HINTS ${MACPORTS_PREFIX}/Library/Frameworks)
    elseif(DEFINED HOMEBREW_PREFIX)
      find_package(OCE 0.17 HINTS ${HOMEBREW_PREFIX}/Cellar/oce/*)
    endif()
  elseif(UNIX)
    set(OCE_DIR "/usr/local/share/cmake/")
  elseif(WIN32)
    set(OCE_DIR "c:/OCE-0.4.0/share/cmake")
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
