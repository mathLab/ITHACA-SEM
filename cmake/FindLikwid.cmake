# FindLikwid
# ----------
# Find likwid, a performance monitoring library.
#
#
# This module will define the following variables:
#
# LIKWID_HEADER_DIR  - Location of likwid headers
# LIKWID_INCLUDE_DIR - Location of likwid includes
# LIKWID_LIBRARY     - Location of likwid dynamic library

IF(NEKTAR_USE_LIKWID)

    FIND_LIBRARY(LIKWID_LIBRARY NAMES likwid PATHS ${LIKWID_DIR}/lib)

    IF(LIKWID_LIBRARY)
        MESSAGE(STATUS "Found likwid library: ${LIKWID_LIBRARY}")
        GET_FILENAME_COMPONENT(LIKWID_LIB_DIR ${LIKWID_LIBRARY} DIRECTORY)
        FIND_PATH(LIKWID_INCLUDE_DIR NAMES likwid.h PATHS ${LIKWID_LIB_DIR}/../include ${LIKWID_DIR}/include)
        IF (LIKWID_INCLUDE_DIR)
            MESSAGE(STATUS "Found likwid header: ${LIKWID_INCLUDE_DIR}")
            ADD_DEFINITIONS(-DLIKWID_PERFMON)

            # INCLUDE_DIRECTORIES(SYSTEM ${LIKWID_INCLUDE_DIR})
            # SET(LIKWID_CONFIG_INCLUDE_DIR ${LIKWID_INCLUDE_DIR})
            # LINK_DIRECTORIES(${LIKWID_LIB_DIR})
            # SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${LIKWID_LIB_DIR}")
            # ADD_CUSTOM_TARGET(likwid ALL)

            MARK_AS_ADVANCED(LIKWID_INCLUDE_DIR)
            MARK_AS_ADVANCED(LIKWID_LIBRARY)
            # MARK_AS_ADVANCED(LIKWID_CONFIG_INCLUDE_DIR)
        ELSE()
            MESSAGE(STATUS "Likwid header not found")
        ENDIF()

    ELSE()
        MESSAGE(STATUS "Likwid library not found")
    ENDIF()

ENDIF(NEKTAR_USE_LIKWID)