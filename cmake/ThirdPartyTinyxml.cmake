# TinyXML
OPTION(NEKTAR_USE_TINYXML_STL "Use STL with TinyXML library." ON)
MARK_AS_ADVANCED(NEKTAR_USE_TINYXML_STL)

OPTION(THIRDPARTY_BUILD_TINYXML 
    "Build TinyXML library from ThirdParty." ON)

IF (THIRDPARTY_BUILD_TINYXML)
    # Tiny XML
    IF (NOT EXISTS ${TPSRC}/tinyxml_2_4_3.tar.bz2)
        FILE(DOWNLOAD ${TPURL}/tinyxml_2_4_3.tar.bz2 
                      ${TPSRC}/tinyxml_2_4_3.tar.bz2)
    ENDIF()
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar jxf ${TPSRC}/tinyxml_2_4_3.tar.bz2
        WORKING_DIRECTORY ${TPSRC}
    )
    SET(TINYXML_INCLUDE_DIR ${TPSRC}/tinyxml)
    SET(TINYXML_BASE ${TPSRC})
    SET(TINYXML_SRC_DIR ${TPSRC}/tinyxml)
ELSE (THIRDPARTY_BUILD_TINYXML)
    INCLUDE (FindTinyXml)
ENDIF (THIRDPARTY_BUILD_TINYXML)
SET(TINYXML_LIB tinyxml)
INCLUDE_DIRECTORIES(SYSTEM ${TINYXML_BASE})



