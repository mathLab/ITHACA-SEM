cmake_minimum_required(VERSION 2.8)

PROJECT(Tinyxml)

SET(TINYXML_SOURCES
    tinyxml.cpp
    tinyxmlparser.cpp
    tinyxmlerror.cpp
    tinystr.cpp
)

SET(TINYXML_HEADERS
    tinystr.h
    tinyxml.h
)

ADD_LIBRARY(tinyxml STATIC ${TINYXML_SOURCES} ${TINYXML_HEADERS})

IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    SET_TARGET_PROPERTIES(tinyxml PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF()

INSTALL(TARGETS tinyxml DESTINATION lib)
INSTALL(FILES ${TINYXML_HEADERS} DESTINATION include)
