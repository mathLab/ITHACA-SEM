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
    FIND_LIBRARY(CCMIO_LIBRARY NAMES "ccmio" "adf" PATHS /usr/local/lib)
    
    IF (CCMIO_LIBRARY)
        MESSAGE(STATUS "Found Ccmio: ${CCMIO_LIBRARY}") 
        MARK_AS_ADVANCED(CCMIO_LIBRARY)
        ADD_DEFINITIONS(-DNEKTAR_USE_CCM)
     ELSE()
        MESSAGE(FATAL_ERROR "Cound not find ccmio library")
     ENDIF()
ENDIF( NEKTAR_USE_CCM )

