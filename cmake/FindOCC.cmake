# Try to find OCE / OCC
# Once done this will define
#
# OCC_FOUND          - system has OCC - OpenCASCADE
# OCC_INCLUDE_DIR    - where the OCC include directory can be found
# OCC_LIBRARY_DIR    - where the OCC library directory can be found
# OCC_LIBRARIES      - Link this to use OCC
# OCC_OCAF_LIBRARIES - Link this to use OCC OCAF framework
#
# Adapted from FreeCAD: http://free-cad.sf.net

set(TEST_ENV $ENV{OCE_ROOT})
if(NOT DEFINED OCE_DIR AND DEFINED TEST_ENV)
  file(GLOB OCE_DIR $ENV{OCE_ROOT}/lib/oce-*)
endif()

set(TEST_ENV $ENV{OCE_DIR})
if(NOT DEFINED OCE_DIR AND DEFINED TEST_ENV)
  set(OCE_DIR $ENV{OCE_DIR})
endif()

SET(OCE_FIND_COMPONENTS
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
    TKXCAF
    TKLCAF
    TKXDESTEP
)

# First try to find OpenCASCADE Community Edition
if(NOT DEFINED OCE_DIR)
  # Check for OSX needs to come first because UNIX evaluates to true on OSX
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if(DEFINED MACPORTS_PREFIX)
      find_package(OCE 0.17 QUIET HINTS ${MACPORTS_PREFIX}/Library/Frameworks)
    elseif(DEFINED HOMEBREW_PREFIX)
      find_package(OCE 0.17 QUIET HINTS ${HOMEBREW_PREFIX}/Cellar/oce/*)
    endif()
  elseif(UNIX)
    set(OCE_DIR "/usr/local/share/cmake/")
    find_package(OCE 0.17 QUIET)
  endif()
endif()

SET(OCC_FOUND FALSE)

if(OCE_FOUND AND OCE_ALL_FOUND)
  message(STATUS "OpenCASCADE Community Edition has been found with all required components.")
  SET(OCC_INCLUDE_DIR ${OCE_INCLUDE_DIRS})
  SET(OCC_FOUND TRUE)
else(OCE_FOUND AND OCE_ALL_FOUND) #look for OpenCASCADE
  message(STATUS "OpenCASCADE Community Edition could not be found or has missing components.")

    FIND_PATH(OCC_INCLUDE_DIR Standard_Version.hxx
      /usr/include/opencascade
      /usr/local/include/opencascade
      /usr/local/opt/opencascade/include
      /opt/opencascade/include
      /opt/opencascade/inc
    )
    FIND_LIBRARY(OCC_LIBRARY TKernel
      /usr/lib
      /usr/local/lib
      /usr/local/opt/opencascade/lib
      /opt/opencascade/lib
      opt/local/lib
    )

  if(OCC_LIBRARY)
    message(STATUS "OpenCASCADE has been found.")
    SET(OCC_FOUND TRUE)
    GET_FILENAME_COMPONENT(OCC_LIBRARY_DIR ${OCC_LIBRARY} PATH)
    IF(NOT OCC_INCLUDE_DIR)
      FIND_PATH(OCC_INCLUDE_DIR Standard_Version.hxx
        ${OCC_LIBRARY_DIR}/../inc
      )
    ENDIF()
  endif(OCC_LIBRARY)
endif(OCE_FOUND AND OCE_ALL_FOUND)

if(OCC_INCLUDE_DIR)
  file(STRINGS ${OCC_INCLUDE_DIR}/Standard_Version.hxx OCC_MAJOR
    REGEX "#define OCC_VERSION_MAJOR.*"
  )
  string(REGEX MATCH "[0-9]+" OCC_MAJOR ${OCC_MAJOR})
  file(STRINGS ${OCC_INCLUDE_DIR}/Standard_Version.hxx OCC_MINOR
    REGEX "#define OCC_VERSION_MINOR.*"
  )
  string(REGEX MATCH "[0-9]+" OCC_MINOR ${OCC_MINOR})
  file(STRINGS ${OCC_INCLUDE_DIR}/Standard_Version.hxx OCC_MAINT
    REGEX "#define OCC_VERSION_MAINTENANCE.*"
  )
  string(REGEX MATCH "[0-9]+" OCC_MAINT ${OCC_MAINT})

  set(OCC_VERSION_STRING "${OCC_MAJOR}.${OCC_MINOR}.${OCC_MAINT}")
endif(OCC_INCLUDE_DIR)

if(OCC_FOUND)
  set(OCC_LIBRARIES ${OCE_FIND_COMPONENTS})
  if(OCC_VERSION_STRING VERSION_LESS 6.8)
    MESSAGE(SEND_ERROR "OCC version too low")
  endif(OCC_VERSION_STRING VERSION_LESS 6.8)
  message(STATUS "-- Found OCE/OpenCASCADE with OCC version: ${OCC_VERSION_STRING}")
endif(OCC_FOUND)
