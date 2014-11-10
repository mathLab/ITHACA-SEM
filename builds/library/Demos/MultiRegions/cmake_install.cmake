# Install script for directory: /tmp/ybao/nektar++/library/Demos/MultiRegions

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

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz1D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz1D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz1D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz1D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz2D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz2D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz3D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz3D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz3DHomo1D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz3DHomo1D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo1D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo1D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz3DHomo2D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Helmholtz3DHomo2D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo2D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo2D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/HDGHelmholtz1D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/HDGHelmholtz1D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz1D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz1D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/HDGHelmholtz2D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/HDGHelmholtz2D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz2D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz2D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/PostProcHDG2D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/PostProcHDG2D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG2D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG2D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/HDGHelmholtz3D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/HDGHelmholtz3D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/PostProcHDG3D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/PostProcHDG3D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG3D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG3D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/HDGHelmholtz3DHomo1D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/HDGHelmholtz3DHomo1D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3DHomo1D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3DHomo1D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Deriv3DHomo1D_SingleMode-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Deriv3DHomo1D_SingleMode"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D_SingleMode-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D_SingleMode"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Deriv3DHomo1D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Deriv3DHomo1D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Deriv3DHomo2D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/Deriv3DHomo2D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo2D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo2D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/SteadyAdvectionDiffusionReaction2D-3.4.0"
    "/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/SteadyAdvectionDiffusionReaction2D"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/SteadyAdvectionDiffusionReaction2D-3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/SteadyAdvectionDiffusionReaction2D"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "demos")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/tmp/ybao/nektar++/builds/library/Demos/MultiRegions/ExtraDemos/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

