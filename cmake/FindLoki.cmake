#
# FindLoki.cmake: Find Loki headers
#

# Try to find system Loki headers
FIND_PATH(LOKI_INCLUDE_DIR loki/Typelist.h)

IF (LOKI_INCLUDE_DIR)
    SET(LOKI_FOUND TRUE)
    SET(LOKI_SYSTEM TRUE)
ENDIF ()

IF (NOT LOKI_FOUND)
    FIND_PATH(LOKI_INCLUDE_DIR loki/Typelist.h
        PATHS ${TPSRC}/loki-0.1.3/include NO_DEFAULT_PATH)

    IF (LOKI_INCLUDE_DIR)
        SET(LOKI_FOUND TRUE)

        # Copy from source to TPDIST
        MESSAGE(STATUS ${LOKI_INCLUDE_DIR})
        FILE(COPY ${LOKI_INCLUDE_DIR}/loki/ DESTINATION ${TPDIST}/include/loki/)
        SET(LOKI_INCLUDE_DIR ${TPDIST}/include CACHE PATH "" FORCE)
    ENDIF()
ENDIF()

IF (LOKI_FOUND)
    MESSAGE(STATUS "Found Loki: ${LOKI_INCLUDE_DIR}")
ELSE()
    MESSAGE(FATAL_ERROR "Could not find Loki")
ENDIF()

MARK_AS_ADVANCED(LOKI_INCLUDE_DIR)
