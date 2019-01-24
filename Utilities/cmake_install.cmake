# Install script for directory: /scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/solvers/IncNavierStokesSolver/Utilities

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist")
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
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CFLStep" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CFLStep")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CFLStep"
         RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/solvers/IncNavierStokesSolver/Utilities/CFLStep")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CFLStep" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CFLStep")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CFLStep"
         OLD_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/ThirdParty/dist/lib:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SolverUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/FieldUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/GlobalMapping:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/MultiRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/Collections:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LocalRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SpatialDomains:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/StdRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LibUtilities:"
         NEW_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CFLStep")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Aliasing" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Aliasing")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Aliasing"
         RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/solvers/IncNavierStokesSolver/Utilities/Aliasing")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Aliasing" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Aliasing")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Aliasing"
         OLD_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/ThirdParty/dist/lib:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SolverUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/FieldUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/GlobalMapping:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/MultiRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/Collections:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LocalRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SpatialDomains:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/StdRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LibUtilities:"
         NEW_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Aliasing")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/NonLinearEnergy" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/NonLinearEnergy")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/NonLinearEnergy"
         RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/solvers/IncNavierStokesSolver/Utilities/NonLinearEnergy")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/NonLinearEnergy" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/NonLinearEnergy")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/NonLinearEnergy"
         OLD_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/ThirdParty/dist/lib:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SolverUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/FieldUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/GlobalMapping:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/MultiRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/Collections:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LocalRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SpatialDomains:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/StdRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LibUtilities:"
         NEW_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/NonLinearEnergy")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Fld2DTo2D5" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Fld2DTo2D5")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Fld2DTo2D5"
         RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/solvers/IncNavierStokesSolver/Utilities/Fld2DTo2D5")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Fld2DTo2D5" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Fld2DTo2D5")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Fld2DTo2D5"
         OLD_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/ThirdParty/dist/lib:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SolverUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/FieldUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/GlobalMapping:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/MultiRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/Collections:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LocalRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SpatialDomains:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/StdRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LibUtilities:"
         NEW_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Fld2DTo2D5")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FldAddFalknerSkanBL" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FldAddFalknerSkanBL")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FldAddFalknerSkanBL"
         RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/solvers/IncNavierStokesSolver/Utilities/FldAddFalknerSkanBL")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FldAddFalknerSkanBL" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FldAddFalknerSkanBL")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FldAddFalknerSkanBL"
         OLD_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/ThirdParty/dist/lib:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SolverUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/FieldUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/GlobalMapping:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/MultiRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/Collections:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LocalRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SpatialDomains:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/StdRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LibUtilities:"
         NEW_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FldAddFalknerSkanBL")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/AddModeTo2DFld" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/AddModeTo2DFld")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/AddModeTo2DFld"
         RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/solvers/IncNavierStokesSolver/Utilities/AddModeTo2DFld")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/AddModeTo2DFld" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/AddModeTo2DFld")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/AddModeTo2DFld"
         OLD_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/ThirdParty/dist/lib:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SolverUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/FieldUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/GlobalMapping:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/MultiRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/Collections:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LocalRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SpatialDomains:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/StdRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LibUtilities:"
         NEW_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/AddModeTo2DFld")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractMeanModeFromHomo1DFld" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractMeanModeFromHomo1DFld")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractMeanModeFromHomo1DFld"
         RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/solvers/IncNavierStokesSolver/Utilities/ExtractMeanModeFromHomo1DFld")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractMeanModeFromHomo1DFld" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractMeanModeFromHomo1DFld")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractMeanModeFromHomo1DFld"
         OLD_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/ThirdParty/dist/lib:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SolverUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/FieldUtils:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/GlobalMapping:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/MultiRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/Collections:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LocalRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/SpatialDomains:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/StdRegions:/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/library/LibUtilities:"
         NEW_RPATH "/scratch/mhess/nektar_cppITHACA/nektar++-4.4.1/dist/lib64/nektar++-4.4.1")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractMeanModeFromHomo1DFld")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "solvers")

