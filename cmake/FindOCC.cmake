SET(TEST_ENV $ENV{OCE_ROOT})
IF(NOT DEFINED OCE_DIR AND DEFINED TEST_ENV)
    FILE(GLOB OCE_DIR $ENV{OCE_ROOT}/lib/oce-*)
ENDIF()

SET(TEST_ENV $ENV{OCE_DIR})
IF(NOT DEFINED OCE_DIR AND DEFINED TEST_ENV)
    SET(OCE_DIR $ENV{OCE_DIR})
ENDIF()

SET(OCE_FIND_COMPONENTS ${OCC_LIB_LIST})

#First try to find OpenCASCADE Community Edition
IF(NOT DEFINED OCE_DIR)
#Check for OSX needs to come first because UNIX evaluates to true on OSX
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        IF(DEFINED MACPORTS_PREFIX)
            FIND_PACKAGE(OCE 0.17 QUIET HINTS ${MACPORTS_PREFIX}/Library/Frameworks)
        ELSEIF(DEFINED HOMEBREW_PREFIX)
            FIND_PACKAGE(OCE 0.17 QUIET HINTS ${HOMEBREW_PREFIX}/Cellar/oce/*)
        ENDIF()
    ENDIF()
ENDIF()

# The version of OCE now being used as the third party version that can be built from
# source as part of the Nektar++ build is 0.18.3. However the FIND_PACKAGE command below
# is looking for a version of OCE that is pre-installed on the system (e.g. from packages)
# The packaged version on Centos 7 is 0.17.1 so this needs to remain at 0.17 to discover
# OCE on Centos 7 (as used in the CI builds). The CMake config scripts installed with the
# packaged version of OCE on Ubuntu will find the version specified here, or a more recent
# version, so this FIND_PACKAGE command successfully results in OCE 0.18.3 installed from
# packages on Ubuntu being discovered.
FIND_PACKAGE(OCE 0.17 QUIET)

SET(OCC_FOUND FALSE)

IF(OCE_FOUND AND OCE_ALL_FOUND)
    MESSAGE(STATUS "OpenCASCADE Community Edition has been found with all required components.")
    SET(OCC_INCLUDE_DIR ${OCE_INCLUDE_DIRS})
    SET(OCC_FOUND TRUE)
ELSE(OCE_FOUND AND OCE_ALL_FOUND) #look for OpenCASCADE
    MESSAGE(STATUS "OpenCASCADE Community Edition could not be found or has missing components.")
    SET(OpenCASCADE_FIND_COMPONENTS ${OCE_FIND_COMPONENTS})
    FIND_PACKAGE(OpenCASCADE)

    IF(OpenCASCADE_FOUND)
        MESSAGE(STATUS "OpenCASCADE has been found with all required components.")
        SET(OCC_INCLUDE_DIR ${OpenCASCADE_INCLUDE_DIR})
        SET(OCC_FOUND TRUE)
    ELSE()
        MESSAGE(STATUS "OpenCASCADE could not be found or has missing components.")
    ENDIF()

ENDIF(OCE_FOUND AND OCE_ALL_FOUND)

IF(OCC_FOUND)
    FILE(STRINGS ${OCC_INCLUDE_DIR}/Standard_Version.hxx OCC_MAJOR
        REGEX "#define OCC_VERSION_MAJOR.*")
    STRING(REGEX MATCH "[0-9]+" OCC_MAJOR ${OCC_MAJOR})
    FILE(STRINGS ${OCC_INCLUDE_DIR}/Standard_Version.hxx OCC_MINOR
        REGEX "#define OCC_VERSION_MINOR.*")
    STRING(REGEX MATCH "[0-9]+" OCC_MINOR ${OCC_MINOR})
    FILE(STRINGS ${OCC_INCLUDE_DIR}/Standard_Version.hxx OCC_MAINT
        REGEX "#define OCC_VERSION_MAINTENANCE.*")
    STRING(REGEX MATCH "[0-9]+" OCC_MAINT ${OCC_MAINT})

    SET(OCC_VERSION_STRING "${OCC_MAJOR}.${OCC_MINOR}.${OCC_MAINT}")

    IF(OCC_VERSION_STRING VERSION_LESS 6.8)
        MESSAGE(STATUS "OCC version too low, will build from source")
        SET(OCC_FOUND FALSE)
    ELSE()
        MESSAGE(STATUS "-- Found OCE/OpenCASCADE with OCC version: ${OCC_VERSION_STRING}")
        SET(OCC_LIBRARIES ${OCE_FIND_COMPONENTS})
    ENDIF()
ENDIF(OCC_FOUND)
