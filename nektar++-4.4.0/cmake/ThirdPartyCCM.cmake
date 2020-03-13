########################################################################
#
# ThirdParty configuration for Nektar++
#
# Star CCM i/o
#
########################################################################
OPTION(NEKTAR_USE_CCM
   "CCM star i/o library is available." OFF)

IF( NEKTAR_USE_CCM )

    set(CCMIO_LIBRARIES
        ccmio
        adf
    )

    FIND_LIBRARY(CCMIO_LIBRARY NAMES "ccmio" PATHS /usr/local/lib ${Nektar++_TP_LIBRARY_DIRS})

    IF( CCMIO_LIBRARY )
        MESSAGE(STATUS "Found Ccmio: ${CCMIO_LIBRARY}")
        MARK_AS_ADVANCED(CCMIO_LIBRARY)
        ADD_DEFINITIONS(-DNEKTAR_USE_CCM)
        FIND_PATH (CCMIO_INCLUDE_DIR ccmio.h)
        GET_FILENAME_COMPONENT(CCMIO_LIBRARY_DIR ${CCMIO_LIBRARY} PATH)
        INCLUDE_DIRECTORIES(NekMesh ${CCMIO_INCLUDE_DIR})
        LINK_DIRECTORIES(${CCMIO_LIBRARY_DIR})
        MESSAGE(STATUS ${CCMIO_LIBRARY_DIR})
     ELSE()
        MESSAGE(FATAL_ERROR "Cound not find ccmio library")
     ENDIF()
ENDIF( NEKTAR_USE_CCM )

