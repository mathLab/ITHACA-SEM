#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
SET(CMAKE_IMPORT_FILE_VERSION 1)

# Compute the installation prefix relative to this file.
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)

# Import target "LibUtilities" for configuration "Release"
SET_PROPERTY(TARGET LibUtilities APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(LibUtilities PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/apps/fftw/3.3.2/lib/libfftw3.so;/apps/boost/1.49.0/lib/libboost_thread.so;/apps/boost/1.49.0/lib/libboost_iostreams.so;/apps/boost/1.49.0/lib/libboost_date_time.so;/apps/boost/1.49.0/lib/libboost_program_options.so;/apps/boost/1.49.0/lib/libboost_filesystem.so;/apps/boost/1.49.0/lib/libboost_system.so;/usr/lib/x86_64-linux-gnu/libz.so;tinyxml;rt;/usr/lib/openmpi/lib/libmpi_cxx.so;/usr/lib/openmpi/lib/libmpi.so;/usr/lib/openmpi/lib/libopen-rte.so;/usr/lib/openmpi/lib/libopen-pal.so;/usr/lib/x86_64-linux-gnu/libdl.so;/usr/lib/x86_64-linux-gnu/libnsl.so;/usr/lib/x86_64-linux-gnu/libutil.so;/usr/lib/x86_64-linux-gnu/libm.so;/usr/lib/x86_64-linux-gnu/libdl.so;/usr/lib/liblapack.so;/usr/lib/libblas.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libLibUtilities.so.3.4.0"
  IMPORTED_SONAME_RELEASE "libLibUtilities.so.3.4.0"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS LibUtilities )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_LibUtilities "${_IMPORT_PREFIX}/lib64/libLibUtilities.so.3.4.0" )

# Import target "LocalRegions" for configuration "Release"
SET_PROPERTY(TARGET LocalRegions APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(LocalRegions PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "SpatialDomains;StdRegions;LibUtilities"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libLocalRegions.so.3.4.0"
  IMPORTED_SONAME_RELEASE "libLocalRegions.so.3.4.0"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS LocalRegions )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_LocalRegions "${_IMPORT_PREFIX}/lib64/libLocalRegions.so.3.4.0" )

# Import target "MultiRegions" for configuration "Release"
SET_PROPERTY(TARGET MultiRegions APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(MultiRegions PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "LocalRegions"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libMultiRegions.so.3.4.0"
  IMPORTED_SONAME_RELEASE "libMultiRegions.so.3.4.0"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS MultiRegions )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_MultiRegions "${_IMPORT_PREFIX}/lib64/libMultiRegions.so.3.4.0" )

# Import target "SpatialDomains" for configuration "Release"
SET_PROPERTY(TARGET SpatialDomains APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(SpatialDomains PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "StdRegions;LibUtilities"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libSpatialDomains.so.3.4.0"
  IMPORTED_SONAME_RELEASE "libSpatialDomains.so.3.4.0"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS SpatialDomains )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_SpatialDomains "${_IMPORT_PREFIX}/lib64/libSpatialDomains.so.3.4.0" )

# Import target "StdRegions" for configuration "Release"
SET_PROPERTY(TARGET StdRegions APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(StdRegions PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "LibUtilities"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libStdRegions.so.3.4.0"
  IMPORTED_SONAME_RELEASE "libStdRegions.so.3.4.0"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS StdRegions )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_StdRegions "${_IMPORT_PREFIX}/lib64/libStdRegions.so.3.4.0" )

# Import target "SolverUtils" for configuration "Release"
SET_PROPERTY(TARGET SolverUtils APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(SolverUtils PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "MultiRegions"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libSolverUtils.so.3.4.0"
  IMPORTED_SONAME_RELEASE "libSolverUtils.so.3.4.0"
  )

LIST(APPEND _IMPORT_CHECK_TARGETS SolverUtils )
LIST(APPEND _IMPORT_CHECK_FILES_FOR_SolverUtils "${_IMPORT_PREFIX}/lib64/libSolverUtils.so.3.4.0" )

# Loop over all imported files and verify that they actually exist
FOREACH(target ${_IMPORT_CHECK_TARGETS} )
  FOREACH(file ${_IMPORT_CHECK_FILES_FOR_${target}} )
    IF(NOT EXISTS "${file}" )
      MESSAGE(FATAL_ERROR "The imported target \"${target}\" references the file
   \"${file}\"
but this file does not exist.  Possible reasons include:
* The file was deleted, renamed, or moved to another location.
* An install or uninstall procedure did not complete successfully.
* The installation package was faulty and contained
   \"${CMAKE_CURRENT_LIST_FILE}\"
but not all the files it references.
")
    ENDIF()
  ENDFOREACH()
  UNSET(_IMPORT_CHECK_FILES_FOR_${target})
ENDFOREACH()
UNSET(_IMPORT_CHECK_TARGETS)

# Cleanup temporary variables.
SET(_IMPORT_PREFIX)

# Commands beyond this point should not need to know the version.
SET(CMAKE_IMPORT_FILE_VERSION)
