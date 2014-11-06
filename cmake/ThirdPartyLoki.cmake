########################################################################
#
# ThirdParty configuration for Nektar++
#
# Loki headers library
#
########################################################################

# Try to find system Loki headers. Hint /opt/local/include for MacPorts
# (although there is no Portfile for Loki currently).
FIND_PATH(LOKI_INCLUDE_DIR loki/Typelist.h PATHS /opt/local/include)

IF (LOKI_INCLUDE_DIR)
    SET(BUILD_LOKI OFF)
ELSE()
    SET(BUILD_LOKI ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_LOKI
    "Download and extract Loki library to ThirdParty." ${BUILD_LOKI})

IF (THIRDPARTY_BUILD_LOKI)
    # Download Loki if it doesn't already exist.
    IF (NOT EXISTS ${TPSRC}/loki-0.1.3.tar.bz2)
        FILE(DOWNLOAD ${TPURL}/loki-0.1.3.tar.bz2 ${TPSRC}/loki-0.1.3.tar.bz2)
    ENDIF()

    # TODO: Check hashes.
    
    # Extract.
    EXECUTE_PROCESS(
        COMMAND ${CMAKE_COMMAND} -E tar jxf ${TPSRC}/loki-0.1.3.tar.bz2
        WORKING_DIRECTORY ${TPSRC}
    )

    # Set LOKI_INCLUDE_DIR.
    FILE(COPY ${TPSRC}/loki-0.1.3/include/loki/ DESTINATION ${TPDIST}/include/loki/)
    SET(LOKI_INCLUDE_DIR ${TPDIST}/include CACHE PATH "" FORCE)

    MESSAGE(STATUS "Build Loki: ${LOKI_INCLUDE_DIR}")
    SET(LOKI_CONFIG_INCLUDE_DIR ${TPINC})
ELSE()
    MESSAGE(STATUS "Found Loki: ${LOKI_INCLUDE_DIR}")
    SET(LOKI_CONFIG_INCLUDE_DIR ${LOKI_INCLUDE_DIR})
ENDIF()

INCLUDE_DIRECTORIES(SYSTEM ${LOKI_INCLUDE_DIR})
MARK_AS_ADVANCED(LOKI_INCLUDE_DIR)
