# TinyXML
OPTION(NEKTAR_USE_TINYXML_STL "Use STL with TinyXML library." ON)
MARK_AS_ADVANCED(NEKTAR_USE_TINYXML_STL)

OPTION(THIRDPARTY_BUILD_TINYXML 
    "Build TinyXML library from ThirdParty." ON)

IF (THIRDPARTY_BUILD_TINYXML)
    #SET(TINYXML_DIR ${TPSRC}/tinyxml)
    #EXTERNALPROJECT_ADD(
    #    tinyxml
    #    PREFIX ${TPSRC}/build
    #    URL ${TPSRC}/tinyxml.tar.bz2
    #    URL_MD5 "aec842139928e65aa7abdff6de0a09ec"
    #    CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX:    PATH=${TPSRC}/   #build/dist ${TPSRC}/build/src/tinyxml
    #)
    #SET(TINYXML_LIB ${TPSRC}/build/dist/lib/libtinyxml.so)
    #SET(TINYXML_BASE ${TPSRC}/build/src)

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
INCLUDE_DIRECTORIES(${TINYXML_BASE})



