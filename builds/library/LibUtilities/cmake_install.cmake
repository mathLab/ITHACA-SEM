# Install script for directory: /tmp/ybao/nektar++/library/LibUtilities

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

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "lib")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY OPTIONAL FILES
    "/tmp/ybao/nektar++/builds/library/LibUtilities/libLibUtilities.so.3.4.0"
    "/tmp/ybao/nektar++/builds/library/LibUtilities/libLibUtilities.so"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libLibUtilities.so.3.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libLibUtilities.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "lib")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "dev")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/ExpressionTemplates" TYPE FILE FILES
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/AssociativeTraits.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/AssociativeTransform.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/BackwardInverseTransform.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/CommutativeTraits.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/CommutativeTransform.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/CreateFromTree.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/ExpressionEvaluator.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/ExpressionTemplates.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/ForwardInverseTransform.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/InverseOperatorTypeTraits.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/InvertNode.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/Node.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/Operators.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/PerformCommutativeTransformIfNeeded.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/PushDownUnaryNodes.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/RemoveAllUnecessaryTemporaries.hpp"
    "/tmp/ybao/nektar++/library/LibUtilities/../ExpressionTemplates/SortAssociativeCommutativeClusters.hpp"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "dev")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "dev")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/LibUtilities" TYPE DIRECTORY FILES "/tmp/ybao/nektar++/library/LibUtilities/./" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\.hpp$")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "dev")

