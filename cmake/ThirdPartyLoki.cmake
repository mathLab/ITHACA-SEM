# Loki
OPTION(THIRDPARTY_BUILD_LOKI
    "Build Loki library from ThirdParty." ON)

IF (THIRDPARTY_BUILD_LOKI)
    IF (NOT EXISTS ${TPSRC}/loki-0.1.3.tar.bz2)
        FILE(DOWNLOAD ${TPURL}/loki-0.1.3.tar.bz2 ${TPSRC}/loki-0.1.3.tar.bz2)
    ENDIF()
    EXECUTE_PROCESS(
        COMMAND ${CMAKE_COMMAND} -E tar jxf ${TPSRC}/loki-0.1.3.tar.bz2
        WORKING_DIRECTORY ${TPSRC}
    )
ENDIF()

INCLUDE (FindLoki)
INCLUDE_DIRECTORIES(SYSTEM ${LOKI_INCLUDE_DIR})
