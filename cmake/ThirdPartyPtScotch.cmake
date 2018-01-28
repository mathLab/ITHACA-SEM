########################################################################
#
# ThirdParty configuration for Nektar++
#
# Scotch partitioner
#
########################################################################

IF (NOT WIN32)
    OPTION(NEKTAR_USE_PTSCOTCH
        "Use PtScotch library for parallel mesh partitioning." OFF)
ENDIF(NOT WIN32)

IF (NEKTAR_USE_PTSCOTCH)
    # First search for system TinyXML installs. Hint /opt/local for MacPorts.
    FIND_LIBRARY(PTSCOTCH_LIBRARY    NAMES ptscotch PATHS ${MACPORTS_PREFIX}/lib)
    FIND_LIBRARY(PTSCOTCHERR_LIBRARY NAMES ptscotcherr PATHS ${MACPORTS_PREFIX}/lib)
    FIND_PATH   (PTSCOTCH_INCLUDE_DIR ptscotch.h PATHS ${MACPORTS_PREFIX}/include /usr/include/scotch)

    MESSAGE(STATUS ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY} ${PTSCOTCH_INCLUDE_DIR})
    
    IF (PTSCOTCH_LIBRARY AND PTSCOTCHERR_LIBRARY AND PTSCOTCH_INCLUDE_DIR)
        SET(BUILD_PTSCOTCH OFF)
    ELSE()
        SET(BUILD_PTSCOTCH ON)
    ENDIF ()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_PTSCOTCH
        "Build Scotch library from ThirdParty" ${BUILD_PTSCOTCH}
        "NEKTAR_USE_PTSCOTCH" OFF)

    IF (THIRDPARTY_BUILD_PTSCOTCH)
        MESSAGE(FATAL_ERROR "build ptscotch not working yet")
    ELSE (THIRDPARTY_BUILD_PTSCOTCH)
        ADD_CUSTOM_TARGET(ptscotch-6.0.4 ALL)
        MESSAGE(STATUS "Found PtScotch: ${PTSCOTCH_LIBRARY}")
        SET(PTSCOTCH_CONFIG_INCLUDE_DIR ${PTSCOTCH_INCLUDE_DIR})
    ENDIF (THIRDPARTY_BUILD_PTSCOTCH)

    INCLUDE_DIRECTORIES(${PTSCOTCH_INCLUDE_DIR})

    MARK_AS_ADVANCED(PTSCOTCH_LIBRARY)
    MARK_AS_ADVANCED(PTSCOTCHERR_LIBRARY)
    MARK_AS_ADVANCED(PTSCOTCH_INCLUDE_DIR)
ENDIF()
