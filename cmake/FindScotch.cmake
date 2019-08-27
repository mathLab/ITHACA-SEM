# FindScotch
# ----------
# Finds the graph partitioning library Scotch.
#
# This module defines the following variables
#   SCOTCH_INCLUDE_DIRS - Location of scotch.h
#   SCOTCH_LIBRARY - Location of main scotch library
#   SCOTCHERR_LIBRARY - Location of scotcherr library
#   SCOTCH_FOUND - TRUE if Scotch found
#   SCOTCH_VERSION - Version of scotch library found
#
# If the component 'ptscotch' is specified
#   PTSCOTCH_LIBRARY - Location of the ptscotch main library
#   PTSCOTCHERR_LIBRARY - Location of the ptscotcherr library
#
# The following environmental variables are used to help find the library:
#   SCOTCH_DIR
#   SCOTCH_INCDIR
#
# Usage:
#   FIND_PACKAGE(Scotch)    - find any version of serial scotch
#   FIND_PACKAGE(Scotch 5)  - find at least version 5 of serial scotch
#   FIND_PACKAGE(Scotch 5 COMPONENTS ptscotch)
#                           - find at least version 5 of PT-scotch

# Determine if we are looking for PT-scotch, or just scotch
SET(PARALLEL OFF)
LIST(FIND Scotch_FIND_COMPONENTS "ptscotch" FIND_PARALLEL)
IF (FIND_PARALLEL GREATER -1)
    SET(PARALLEL ON)
ENDIF()

MESSAGE(STATUS "Searching for Scotch:")
SET(TEST_SCOTCH_DIR $ENV{SCOTCH_DIR})
SET(TEST_SCOTCH_HOME $ENV{SCOTCH_HOME})
SET(TEST_SCOTCH_INCLUDE_DIR $ENV{SCOTCH_INCDIR})

SET(SCOTCH_HEADERS_DIRS "SCOTCH_HEADERS_DIR-NOTFOUND")
IF(TEST_SCOTCH_INCLUDE_DIR)
    FIND_PATH(SCOTCH_HEADERS_DIRS NAMES scotch.h
                HINTS ${TEST_SCOTCH_INCLUDE_DIR})
ELSEIF(TEST_SCOTCH_DIR)
    FIND_PATH(SCOTCH_HEADERS_DIRS NAMES scotch.h
                HINTS ${TEST_SCOTCH_DIR}
                PATH_SUFFIXES "include" "include/scotch")
ELSEIF(TEST_SCOTCH_HOME)
    FIND_PATH(SCOTCH_HEADERS_DIRS NAMES scotch.h
                HINTS ${TEST_SCOTCH_HOME}
                PATH_SUFFIXES "include" "include/scotch")
ELSE()
    FIND_PATH(SCOTCH_HEADERS_DIRS NAMES scotch.h
                HINTS ${MACPORTS_PREFIX}/include
                PATH_SUFFIXES "scotch")
ENDIF()
MARK_AS_ADVANCED(SCOTCH_HEADERS_DIRS)

IF (SCOTCH_HEADERS_DIRS)
    SET(SCOTCH_INCLUDE_DIR ${SCOTCH_HEADERS_DIRS})

    TRY_RUN(
        RUN_RESULT COMPILE_RESULT
        ${CMAKE_CURRENT_BINARY_DIR}/
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/scripts/get-scotch-version.c
        CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${SCOTCH_INCLUDE_DIR}
        COMPILE_OUTPUT_VARIABLE COMPILER_OUTPUT
        RUN_OUTPUT_VARIABLE SCOTCH_VERSION)

    IF("${COMPILE_RESULT}" AND ("${RUN_RESULT}" EQUAL 0))
        STRING(STRIP "${SCOTCH_VERSION}" SCOTCH_VERSION)
    ELSE()
        MESSAGE(STATUS "${COMPILER_OUTPUT}")
        MESSAGE(STATUS "${SCOTCH_VERSION}")
        MESSAGE(ERROR "ERROR checking for Scotch version.")
    ENDIF()
ELSE ()
    SET(SCOTCH_INCLUDE_DIR "SCOTCH_INCLUDE_DIRS-NOTFOUND")
    IF (NOT SCOTCH_FIND_QUIETLY)
        MESSAGE(STATUS "Looking for Scotch -- scotch.h not found")
    ENDIF ()
ENDIF ()
LIST(REMOVE_DUPLICATES SCOTCH_INCLUDE_DIR)

# Search for the library also in the ../lib directory of scotch.h
GET_FILENAME_COMPONENT(SEARCH_PATHS ${SCOTCH_INCLUDE_DIR} DIRECTORY)
FIND_LIBRARY(SCOTCH_LIBRARY    NAMES scotch PATHS ${SEARCH_PATHS}
    PATH_SUFFIXES lib)
FIND_LIBRARY(SCOTCHERR_LIBRARY NAMES scotcherr PATHS ${SEARCH_PATHS}
    PATH_SUFFIXES lib)
GET_FILENAME_COMPONENT(SCOTCH_LIBRARY_DIR ${SCOTCH_LIBRARY} PATH)

IF (SCOTCH_LIBRARY AND SCOTCHERR_LIBRARY AND SCOTCH_INCLUDE_DIR)
    SET(Scotch_scotch_FOUND TRUE)

    IF (PARALLEL)
        # Start from clean slate
        UNSET(PTSCOTCH_LIBRARY CACHE)

        IF (DEFINED ENV{LD_LIBRARY_PATH})
            STRING(REPLACE ":" ";" SEARCH_LIB_PATH $ENV{LD_LIBRARY_PATH})
        ENDIF()
        FIND_LIBRARY(PTSCOTCH_LIBRARY NAMES ptscotch PATHS ${SCOTCH_LIBRARY_DIR}
            ${SEARCH_LIB_PATH} NO_DEFAULT_PATH)
        FIND_LIBRARY(PTSCOTCHERR_LIBRARY NAMES ptscotcherr PATHS
            ${SCOTCH_LIBRARY_DIR} ${SEARCH_LIB_PATH}
            NO_DEFAULT_PATH)

        IF (PTSCOTCH_LIBRARY AND PTSCOTCHERR_LIBRARY)
            SET(Scotch_ptscotch_FOUND TRUE)

            # Finally, test if serial library needs to be linked as well.
            TRY_COMPILE( COMPILE_RESULT
                ${CMAKE_CURRENT_BINARY_DIR}/
                ${CMAKE_CURRENT_SOURCE_DIR}/cmake/scripts/check-ptscotch-link.c
                CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${SCOTCH_INCLUDE_DIR}
                    -DLINK_DIRECTORIES=${SCOTCH_LIBRARY_DIR}
                LINK_LIBRARIES "ptscotch" "ptscotcherr"
                OUTPUT_VARIABLE COMPILER_OUTPUT)

            IF (COMPILE_RESULT)
                MESSAGE(STATUS "Guessing that serial Scotch does not need to be linked")
                SET(IS_SERIAL_NEEDED FALSE)
            ELSE()
                MESSAGE(STATUS "Guessing that serial Scotch needs to be linked as well")
                SET(IS_SERIAL_NEEDED TRUE)
            ENDIF()

            IF (IS_SERIAL_NEEDED)
                GET_FILENAME_COMPONENT(PTSCOTCH_BASE_DIR ${PTSCOTCH_LIBRARY} DIRECTORY)
                FIND_LIBRARY(PTSCOTCH_LIBRARY_SERIAL NAMES scotch
                    PATHS ${PTSCOTCH_BASE_DIR} PTSCOTCH_LIBRARY_SERIAL NO_DEFAULT_PATH)
                MARK_AS_ADVANCED(PTSCOTCH_LIBRARY_SERIAL)
                IF (PTSCOTCH_LIBRARY_SERIAL)
                    SET(PTSCOTCH_LIBRARY ${PTSCOTCH_LIBRARY} ${PTSCOTCH_LIBRARY_SERIAL}
                        CACHE FILEPATH "PtScotch library" FORCE)
                ENDIF()
            ENDIF()
        ENDIF()
    ENDIF()
ENDIF()

SET(SCOTCH_FOUND TRUE)
IF (Scotch_scotch_FOUND)
    MESSAGE(STATUS "-- Found Scotch version ${SCOTCH_VERSION}: ${SCOTCH_LIBRARY}")
    IF (PARALLEL)
        IF (Scotch_ptscotch_FOUND)
            MESSAGE(STATUS "-- Found PT-Scotch version ${SCOTCH_VERSION}:${PTSCOTCH_LIBRARY}")
        ELSE()
            SET(SCOTCH_FOUND FALSE)
            MESSAGE(STATUS "-- Could not find PT-Scotch library")
        ENDIF()
    ENDIF()
    IF (Scotch_FIND_VERSION)
        IF (${Scotch_FIND_VERSION} VERSION_GREATER ${SCOTCH_VERSION})
            SET(SCOTCH_FOUND FALSE)
            MESSAGE(STATUS "-- This scotch is too old")
        ENDIF()
    ENDIF()
ELSE()
    SET(SCOTCH_FOUND FALSE)
    MESSAGE(STATUS "-- Could not find Scotch library")
ENDIF()
