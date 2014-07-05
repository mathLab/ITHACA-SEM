# TinyXML
OPTION(NEKTAR_USE_TINYXML_STL "Use STL with TinyXML library." ON)
MARK_AS_ADVANCED(NEKTAR_USE_TINYXML_STL)

OPTION(THIRDPARTY_BUILD_TINYXML 
    "Build TinyXML library from ThirdParty." ON)

IF (THIRDPARTY_BUILD_TINYXML)
    INCLUDE(ExternalProject)
    EXTERNALPROJECT_ADD(
        tinyxml-2.4.3
        PREFIX ${TPSRC}
        URL ${TPURL}/tinyxml_2_4_3-1.tar.bz2
        URL_MD5 "76e7fa1264520ec0bbc36a0d5117806d"
        DOWNLOAD_DIR ${TPSRC}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPSRC}/dist
            -DCMAKE_CXX_FLAGS:STRING=-DTIXML_USE_STL
            ${TPSRC}/src/tinyxml-2.4.3
        INSTALL_COMMAND $(MAKE) install
            COMMAND ${CMAKE_COMMAND} -E copy 
                ${TPSRC}/src/tinyxml-2.4.3/tinyxml.h
                ${TPSRC}/dist/include/tinyxml.h
    )
    SET(TINYXML_LIB tinyxml CACHE FILEPATH 
        "Tinyxml library" FORCE)
    MARK_AS_ADVANCED(TINYXML_LIB)
    LINK_DIRECTORIES(${TPSRC}/dist/lib)
    INCLUDE_DIRECTORIES(${TPSRC}/dist/include)
ELSE (THIRDPARTY_BUILD_TINYXML)
    INCLUDE (FindTinyXml)
ENDIF (THIRDPARTY_BUILD_TINYXML)




