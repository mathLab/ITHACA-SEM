# Install script for directory: /tmp/ybao/nektar++

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/tmp/ybao/nektar++/builds/dist")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/Nektar++Libraries.cmake")
    FILE(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/Nektar++Libraries.cmake"
         "/tmp/ybao/nektar++/builds/CMakeFiles/Export/lib64/Nektar++Libraries.cmake")
    IF(EXPORT_FILE_CHANGED)
      FILE(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/Nektar++Libraries-*.cmake")
      IF(OLD_CONFIG_FILES)
        MESSAGE(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/Nektar++Libraries.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        FILE(REMOVE ${OLD_CONFIG_FILES})
      ENDIF(OLD_CONFIG_FILES)
    ENDIF(EXPORT_FILE_CHANGED)
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE FILE FILES "/tmp/ybao/nektar++/builds/CMakeFiles/Export/lib64/Nektar++Libraries.cmake")
  IF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE FILE FILES "/tmp/ybao/nektar++/builds/CMakeFiles/Export/lib64/Nektar++Libraries-release.cmake")
  ENDIF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/tmp/ybao/nektar++/builds/dist/Nektar++Config.cmake;/tmp/ybao/nektar++/builds/dist/FindLoki.cmake;/tmp/ybao/nektar++/builds/dist/FindAccelerateFramework.cmake;/tmp/ybao/nektar++/builds/dist/FindCHUDFramework.cmake;/tmp/ybao/nektar++/builds/dist/FindACML.cmake;/tmp/ybao/nektar++/builds/dist/FindArpack.cmake;/tmp/ybao/nektar++/builds/dist/FindOpenBlas.cmake;/tmp/ybao/nektar++/builds/dist/FindNativeBlasLapack.cmake;/tmp/ybao/nektar++/builds/dist/FindMKL.cmake;/tmp/ybao/nektar++/builds/dist/FindMetis.cmake;/tmp/ybao/nektar++/builds/dist/FindScotch.cmake;/tmp/ybao/nektar++/builds/dist/FindFFTW.cmake;/tmp/ybao/nektar++/builds/dist/FindWin32Lapack.cmake;/tmp/ybao/nektar++/builds/dist/FindTinyXml.cmake;/tmp/ybao/nektar++/builds/dist/FindGSMPI.cmake;/tmp/ybao/nektar++/builds/dist/FindXXT.cmake;/tmp/ybao/nektar++/builds/dist/FindSMV.cmake")
FILE(INSTALL DESTINATION "/tmp/ybao/nektar++/builds/dist" TYPE FILE FILES
    "/tmp/ybao/nektar++/builds/Nektar++Config.cmake"
    "/tmp/ybao/nektar++/cmake/FindLoki.cmake"
    "/tmp/ybao/nektar++/cmake/FindAccelerateFramework.cmake"
    "/tmp/ybao/nektar++/cmake/FindCHUDFramework.cmake"
    "/tmp/ybao/nektar++/cmake/FindACML.cmake"
    "/tmp/ybao/nektar++/cmake/FindArpack.cmake"
    "/tmp/ybao/nektar++/cmake/FindOpenBlas.cmake"
    "/tmp/ybao/nektar++/cmake/FindNativeBlasLapack.cmake"
    "/tmp/ybao/nektar++/cmake/FindMKL.cmake"
    "/tmp/ybao/nektar++/cmake/FindMetis.cmake"
    "/tmp/ybao/nektar++/cmake/FindScotch.cmake"
    "/tmp/ybao/nektar++/cmake/FindFFTW.cmake"
    "/tmp/ybao/nektar++/cmake/FindWin32Lapack.cmake"
    "/tmp/ybao/nektar++/cmake/FindTinyXml.cmake"
    "/tmp/ybao/nektar++/cmake/FindGSMPI.cmake"
    "/tmp/ybao/nektar++/cmake/FindXXT.cmake"
    "/tmp/ybao/nektar++/cmake/FindSMV.cmake"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/tmp/ybao/nektar++/builds/library/cmake_install.cmake")
  INCLUDE("/tmp/ybao/nektar++/builds/solvers/cmake_install.cmake")
  INCLUDE("/tmp/ybao/nektar++/builds/utilities/cmake_install.cmake")
  INCLUDE("/tmp/ybao/nektar++/builds/tests/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/tmp/ybao/nektar++/builds/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/tmp/ybao/nektar++/builds/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
